import computeShaderWGSL from './shaders/compute.wgsl?raw';
import postfxShaderWGSL from './shaders/postfx.wgsl?raw';
import galaxyShaderWGSL from './shaders/galaxy.wgsl?raw';
import { AppState, CameraMode } from './state';
import orbitsShaderWGSL from './shaders/orbits.wgsl?raw';
import { themes } from './theme';
import { G } from './shared/constants';
import { spectralResponses } from './spectral';
import type { Authority } from './authority/authority';
import type { Body, Ship } from './shared/types';
import { InputManager } from './input';
import type { Camera } from './camera';
import { CameraManager } from './camera/manager';
import { HUDManager } from './hud';
import { Galaxy } from './galaxy';
import type { Star } from './galaxy';
import type { Vec3 } from './shared/types';
import { mat4 } from 'gl-matrix';

// Determine parent-child relationships by strongest gravitational influence
function buildSystemHierarchy(bodies: Body[]): Map<string, string | null> {
  const parentMap = new Map<string, string | null>();
  for (let i = 0; i < bodies.length; i++) {
    const body = bodies[i];
    let bestParentId: string | null = null;
    let maxForce = -Infinity;
    for (let j = 0; j < bodies.length; j++) {
      if (i === j) continue;
      const potentialParent = bodies[j];
      const dx = potentialParent.position[0] - body.position[0];
      const dy = potentialParent.position[1] - body.position[1];
      const dz = potentialParent.position[2] - body.position[2];
      const distSq = dx * dx + dy * dy + dz * dz;
      if (distSq <= 0) continue;
      // G is a constant factor and can be omitted for comparison; included for completeness
      const force = potentialParent.mass / distSq; // (G * m_parent * m_body) / r^2 -> compare by m_parent / r^2
      if (force > maxForce) {
        maxForce = force;
        bestParentId = potentialParent.id;
      }
    }
    parentMap.set(body.id, bestParentId);
  }
  return parentMap;
}

const VISUAL_SETTINGS = {
  systemViewSize: 20.0,
} as const;

export class Renderer {

  private canvas: HTMLCanvasElement;
  private hudCanvas?: HTMLCanvasElement;
  private device!: GPUDevice;
  private context!: GPUCanvasContext;
  private presentationFormat!: GPUTextureFormat;
  private textureSize!: { width: number, height: number };
  private computeShader!: { pipeline: GPUComputePipeline, bindGroup: GPUBindGroup, texture: GPUTexture };
  private renderPipeline!: GPURenderPipeline;
  private presentPipeline!: GPURenderPipeline;
  private sampler!: GPUSampler;
  private themeUniformBuffer!: GPUBuffer;
  public authority: Authority;
  private spheresBuffer!: GPUBuffer;
  private cameraUniformBuffer!: GPUBuffer;
  private cameraManager!: CameraManager;
  private postFxTextureA!: GPUTexture;
  private postFxTextureB!: GPUTexture;
  private frameCount = 0;
  private currentThemeName: string = 'amber';
  private currentResponseName: string = 'Visible (Y)';
  private lastDeltaTime = 1/60;
  private state: AppState;
  private galaxyPipeline!: GPURenderPipeline;
  private starBuffer!: GPUBuffer;
  private galaxyCameraUniformBuffer!: GPUBuffer;
  private galaxy!: Galaxy;
  private galaxyBindGroup!: GPUBindGroup;
  private galaxyDepthTexture!: GPUTexture;
  private orbitsTexture!: GPUTexture;
  // Orbits overlay
  private orbitsPipeline!: GPURenderPipeline;
  private orbitsUniformBuffer!: GPUBuffer;
  private orbitsBindGroup!: GPUBindGroup;
  private targetInfoUniform!: GPUBuffer;
  private orbitBuffers: GPUBuffer[] = [];
  private orbitCPUData: Float32Array[] = [];
  private orbitCounts: number[] = [];
  private readonly ORBIT_MAX_POINTS = 2048;
  private readonly ORBIT_TAIL_HIDE = 0;
  private _lastBodiesForMask: Body[] = [];
  private readonly ORBIT_SAMPLES = 256;
  private lastSystemScale: number = 1.0;
  private lastKnownBodyCount = 0;
  private input: InputManager;
  private hud!: HUDManager;

  constructor(canvas: HTMLCanvasElement, authority: Authority, state: AppState, input: InputManager) {
    this.canvas = canvas;
    this.authority = authority;
    this.cameraManager = new CameraManager(this.canvas);
    this.state = state;
    this.input = input;
    // HUD setup
    const hudCanvas = document.getElementById('hud-canvas') as HTMLCanvasElement | null;
    if (hudCanvas) {
      this.hudCanvas = hudCanvas;
      this.hud = new HUDManager(hudCanvas);
    } else {
      console.warn('HUD canvas not found');
    }
  }

  public getCamera(): Camera {
    return this.cameraManager.getCamera();
  }

  private get camera() {
    return this.cameraManager.getCamera();
  }

  public setTheme(themeName: string, responseName: string, deltaTime?: number): void {
    this.currentThemeName = themeName;
    this.currentResponseName = responseName;
    if (deltaTime !== undefined) this.lastDeltaTime = deltaTime;

    const theme = themes[this.currentThemeName as keyof typeof themes];
    const response = spectralResponses[this.currentResponseName];

    if (!theme || !response) {
      console.error("Invalid theme or response name");
      return;
    }

    // 80-byte buffer, write first 68 bytes used; pad the rest implicitly
    const themeData = new Float32Array(20);
    themeData.set(theme.bg, 0);
    themeData.set(theme.fg, 4);
    themeData.set(theme.accent, 8);
    themeData.set(response, 12);
    themeData[16] = this.lastDeltaTime;
    this.device.queue.writeBuffer(this.themeUniformBuffer, 0, themeData);
  }

  private async initWebGPU() {
    if (!navigator.gpu) {
      throw new Error("WebGPU not supported on this browser.");
    }

    const adapter = await navigator.gpu.requestAdapter();
    if (!adapter) {
      throw new Error("No appropriate GPUAdapter found.");
    }

    const device = await adapter.requestDevice();
    const context = this.canvas.getContext("webgpu") as any;
    const presentationFormat = navigator.gpu.getPreferredCanvasFormat();

    context.configure({
      device,
      format: presentationFormat,
    });

    return { device, context, presentationFormat };
  }

  private createComputeShader(textureSize: { width: number, height: number }, numBodies: number, spheresBuffer: GPUBuffer) {
    const module = this.device.createShaderModule({
      label: "Compute shader",
      code: `const NUM_SPHERES: u32 = ${numBodies};\n` + computeShaderWGSL
    });

    const pipeline = this.device.createComputePipeline({
      label: "Compute pipeline",
      layout: "auto",
      compute: {
        module,
        entryPoint: "main"
      }
    });

    const texture = this.device.createTexture({
      size: textureSize,
      format: 'rgba16float',
      usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING
    });

    const bindGroup = this.device.createBindGroup({
      label: "Compute bind group",
      layout: pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: spheresBuffer } },
        { binding: 1, resource: texture.createView() },
        { binding: 2, resource: { buffer: this.cameraUniformBuffer } }
      ]
    });

    return { pipeline, bindGroup, texture };
  }

  private _reinitializeComputeResources(numBodies: number) {
    const floatsPerSphere = 16;
    if (this.spheresBuffer) {
      try { this.spheresBuffer.destroy(); } catch {}
    }
    this.spheresBuffer = this.device.createBuffer({
      size: Math.max(1, numBodies * floatsPerSphere * 4),
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.computeShader = this.createComputeShader(this.textureSize, numBodies, this.spheresBuffer);
    this.lastKnownBodyCount = numBodies;
  }

  private createPostFXPipeline(format: GPUTextureFormat) {
    const module = this.device.createShaderModule({
      label: "PostFX shader",
      code: postfxShaderWGSL
    });

    return this.device.createRenderPipeline({
      layout: "auto",
      vertex: {
        module,
        entryPoint: "vertexMain"
      },
      fragment: {
        module,
        entryPoint: "fragmentMain",
        targets: [{ format }]
      },
      primitive: {
        topology: "triangle-list"
      }
    });
  }

  private createPresentPipeline(format: GPUTextureFormat) {
    const module = this.device.createShaderModule({
      label: "Present shader",
      code: postfxShaderWGSL
    });

    return this.device.createRenderPipeline({
      layout: "auto",
      vertex: {
        module,
        entryPoint: "vertexMain"
      },
      fragment: {
        module,
        entryPoint: "presentFragment",
        targets: [{ format }]
      },
      primitive: {
        topology: "triangle-list"
      }
    });
  }

  // Clear a renderable texture to a solid color (defaults to opaque black)
  private _clearTexture(texture: GPUTexture, color: { r: number, g: number, b: number, a: number } = { r: 0, g: 0, b: 0, a: 1 }): void {
    const encoder = this.device.createCommandEncoder();
    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: texture.createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: color,
      }]
    });
    pass.end();
    this.device.queue.submit([encoder.finish()]);
  }

  public async start() {
    const { device, context, presentationFormat } = await this.initWebGPU();
    this.device = device;
    this.context = context;
    this.presentationFormat = presentationFormat;
    // Galaxy resources
    this.galaxy = new Galaxy(1337, 20000, 200);
    const starData = this._serializeStars(this.galaxy.stars);
    this.starBuffer = this.device.createBuffer({
      size: starData.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(this.starBuffer, 0, starData.buffer as ArrayBuffer, 0, starData.byteLength);

    // 128 bytes for alignment (64 for viewProjection, 16 for right vec3 padded, 16 for up vec3 padded)
    this.galaxyCameraUniformBuffer = this.device.createBuffer({
      size: 128,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });


    this.cameraUniformBuffer = this.device.createBuffer({
      size: 256,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    // Defer initial camera uniform write until system scale is known

    const width = this.canvas.getBoundingClientRect().width;
    const aspectRation = 16.0 / 9.0;
    const height = Math.floor(width / aspectRation);
    this.canvas.width = Math.max(1, Math.min(width, this.device.limits.maxTextureDimension2D));
    this.canvas.height = Math.max(1, Math.min(height, this.device.limits.maxTextureDimension2D));
    this.textureSize = { width: this.canvas.width, height: this.canvas.height };
    // Sync HUD canvas size
    if (this.hudCanvas) {
      this.hudCanvas.width = this.canvas.width;
      this.hudCanvas.height = this.canvas.height;
    }

    const systemState = await this.authority.query();
    const numBodies = systemState.bodies.length;
    this._reinitializeComputeResources(numBodies);
    const initSerialized = this._serializeSystemState(systemState.bodies, this.camera.focusBodyId);
    const sphereData = initSerialized.sphereData;
    this.device.queue.writeBuffer(this.spheresBuffer, 0, sphereData.buffer as ArrayBuffer, 0, sphereData.byteLength);
    // Initialize camera uniforms: world is camera-relative, so camera is at origin
    {
      const cam = this.camera;
      const cameraData = new Float32Array(16);
      cameraData.set([0, 0, 0], 0);
      const lookAtRelative: Vec3 = [
        cam.look_at[0] - cam.eye[0],
        cam.look_at[1] - cam.eye[1],
        cam.look_at[2] - cam.eye[2],
      ];
      cameraData.set(lookAtRelative, 4);
      cameraData.set(cam.up, 8);
      cameraData[12] = 0.0;
      this.device.queue.writeBuffer(this.cameraUniformBuffer, 0, cameraData);
    }
    this.renderPipeline = this.createPostFXPipeline(this.presentationFormat);
    this.presentPipeline = this.createPresentPipeline(this.presentationFormat);

    // Orbits pipeline & buffers
    this._createOrbitsPipeline();
    // 80 bytes for alignment (64 for viewProjection + 16 padding, store color in first 3 of next vec4)
    this.orbitsUniformBuffer = this.device.createBuffer({
      size: 80,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    // Target info uniform: struct size is 272 bytes (matches WGSL OrbitMask)
    this.targetInfoUniform = this.device.createBuffer({
      size: 272,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.orbitsBindGroup = this.device.createBindGroup({
      layout: this.orbitsPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: this.orbitsUniformBuffer } },
      ],
    });
    // One vertex buffer per body
    this.orbitBuffers = systemState.bodies.map(() => this.device.createBuffer({
      size: this.ORBIT_MAX_POINTS * 3 * 4,
      usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
    }));
    this.orbitCPUData = systemState.bodies.map(() => new Float32Array(this.ORBIT_MAX_POINTS * 3));
    this.orbitCounts = systemState.bodies.map(() => 0);
    this._lastBodiesForMask = systemState.bodies;

    // Orbits render target
    this.orbitsTexture = this.device.createTexture({
      size: this.textureSize,
      format: this.presentationFormat,
      usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST,
    });
    // Ensure orbits texture starts cleared so sampling adds no unintended brightness
    this._clearTexture(this.orbitsTexture);

    // Galaxy pipeline and bind group
    this._createGalaxyPipeline();
    this.galaxyBindGroup = this.device.createBindGroup({
      layout: this.galaxyPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: this.galaxyCameraUniformBuffer } },
        { binding: 1, resource: { buffer: this.starBuffer } },
      ],
    });

    this.sampler = this.device.createSampler({
      magFilter: "linear",
      minFilter: "linear"
    });

    // Theme uniform buffer (4 x vec3 padded + f32 deltaTime). Use 80 bytes for alignment simplicity
    this.themeUniformBuffer = this.device.createBuffer({
      size: 80,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    // Upload default theme and response
    this.setTheme(this.currentThemeName, this.currentResponseName, this.lastDeltaTime);

    // Create ping-pong textures for persistence (match presentation format for compatibility)
    this.postFxTextureA = this.device.createTexture({
      size: this.textureSize,
      format: this.presentationFormat,
      usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC,
    });
    this.postFxTextureB = this.device.createTexture({
      size: this.textureSize,
      format: this.presentationFormat,
      usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC,
    });
    // Clear both ping-pong textures so the first frame's prev-frame sampling is black
    this._clearTexture(this.postFxTextureA);
    this._clearTexture(this.postFxTextureB);

    // Depth texture for galaxy raster path
    this.galaxyDepthTexture = this.device.createTexture({
      size: this.textureSize,
      format: 'depth24plus',
      usage: GPUTextureUsage.RENDER_ATTACHMENT,
    });

    const observer = new ResizeObserver(async entries => {
      for (const entry of entries) {
        const renderTime = performance.now();
        const canvas = entry.target as HTMLCanvasElement;
        canvas.width = Math.max(1, Math.min(width, this.device.limits.maxTextureDimension2D));
        canvas.height = Math.max(1, Math.min(height, this.device.limits.maxTextureDimension2D));
        // Resize HUD canvas to match
        if (this.hudCanvas) {
          this.hudCanvas.width = canvas.width;
          this.hudCanvas.height = canvas.height;
        }

        this.context.configure({
          device: this.device,
          format: this.presentationFormat
        });

        this.textureSize = { width: canvas.width, height: canvas.height };
        // Recreate compute pipeline with up-to-date body count
        const currentState = await this.authority.query();
        const correctNumBodies = currentState.bodies.length;
        this.computeShader = this.createComputeShader(this.textureSize, correctNumBodies, this.spheresBuffer);
        // Recreate ping-pong textures on resize
        this.postFxTextureA = this.device.createTexture({
          size: this.textureSize,
          format: this.presentationFormat,
          usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC,
        });
        this.postFxTextureB = this.device.createTexture({
          size: this.textureSize,
          format: this.presentationFormat,
          usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC,
        });
        // Clear ping-pong targets after resize to avoid sampling uninitialized data
        this._clearTexture(this.postFxTextureA);
        this._clearTexture(this.postFxTextureB);
        // Recreate galaxy depth texture on resize
        this.galaxyDepthTexture = this.device.createTexture({
          size: this.textureSize,
          format: 'depth24plus',
          usage: GPUTextureUsage.RENDER_ATTACHMENT,
        });
        // Recreate orbits render target on resize
        this.orbitsTexture = this.device.createTexture({
          size: this.textureSize,
          format: this.presentationFormat,
          usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST,
        });
        this._clearTexture(this.orbitsTexture);
        this._render();
        console.log("Rerender time:", performance.now() - renderTime);
      }
    });
    observer.observe(this.canvas);

    this._update();
  }

  private _serializeSystemState(bodies: Body[], focusBodyId: string | null, fixedScale?: number): { sphereData: Float32Array, focusBodyRenderedRadius: number, systemScale: number } {
    const floatsPerSphere = 16; // 64 bytes per sphere to match WGSL layout
    const sphereData = new Float32Array(bodies.length * floatsPerSphere);
    const sphereDataU32 = new Uint32Array(sphereData.buffer);
    const eye: Vec3 = this.camera.eye;

    let focusBodyRenderedRadius = 0;

    // Determine hierarchy and primary center
    const hierarchy = buildSystemHierarchy(bodies);
    const systemCenters = bodies.filter(b => hierarchy.get(b.id) === null);
    const primaryCenter = systemCenters.length > 0
      ? systemCenters.reduce((a, b) => a.mass > b.mass ? a : b)
      : bodies.reduce((a, b) => a.mass > b.mass ? a : b);

    // Compute system scale normalization (will be overridden if fixedScale provided)
    const systemRadius = bodies.reduce((maxDist, b) => {
      if (b.id === primaryCenter.id) return maxDist;
      const dx = b.position[0] - primaryCenter.position[0];
      const dy = b.position[1] - primaryCenter.position[1];
      const dz = b.position[2] - primaryCenter.position[2];
      const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
      return Math.max(maxDist, dist);
    }, 0);
    let systemScale = systemRadius > 0 ? VISUAL_SETTINGS.systemViewSize / systemRadius : 1.0;
    if (fixedScale !== undefined && fixedScale > 0) {
      systemScale = fixedScale;
    }

    // Determine current viewing context

    const renderScale = systemScale;

    bodies.forEach((body, i) => {
      const f_base = i * floatsPerSphere; // float index base
      const u_base = f_base;              // u32 index base

      const renderedRadius = body.radius * renderScale;
      
      if (body.id === focusBodyId) {
        focusBodyRenderedRadius = renderedRadius;
      }
      if (body.id === focusBodyId) {
        focusBodyRenderedRadius = renderedRadius;
      }

      // Sphere.center (vec3f) in camera-relative coordinates -> slots 0, 1, 2
      sphereData[f_base + 0] = (body.position[0] * renderScale) - eye[0];
      sphereData[f_base + 1] = (body.position[1] * renderScale) - eye[1];
      sphereData[f_base + 2] = (body.position[2] * renderScale) - eye[2];

      // Sphere.radius (f32) -> slot 3
      sphereData[f_base + 3] = renderedRadius;

      // Material starts at slot 4 (offset 16 bytes)
      // material.albedo (vec3f) -> slots 4, 5, 6
      sphereData[f_base + 4] = body.albedo[0];
      sphereData[f_base + 5] = body.albedo[1];
      sphereData[f_base + 6] = body.albedo[2];

      // padding at slot 7 to align next vec3 to 32 bytes

      // material.emissive (vec3f) -> slots 8, 9, 10
      const emissive = body.emissive ?? [0, 0, 0];
      sphereData[f_base + 8] = emissive[0];
      sphereData[f_base + 9] = emissive[1];
      sphereData[f_base + 10] = emissive[2];

      // material.fuzziness (f32) -> slot 11 (offset 44 bytes)
      sphereData[f_base + 11] = 0.0;

      // material.refraction_index (f32) -> slot 12 (offset 48 bytes)
      sphereData[f_base + 12] = 1.0;

      // material.mat_type (u32) -> slot 13 (offset 52 bytes)
      sphereDataU32[u_base + 13] = 0; // Lambertian

      // slots 14, 15 are padding to reach 64 bytes
    });
    return { sphereData, focusBodyRenderedRadius, systemScale };
  }

  private async _update() {
    const deltaTime = 1 / 60.0;
    const inputState = { deltaX: this.input.deltaX, deltaY: this.input.deltaY, keys: this.input.keys };
    await this.authority.tick(deltaTime, inputState);
    const systemState = await this.authority.query();
    if (systemState.bodies.length !== this.lastKnownBodyCount) {
      this._reinitializeComputeResources(systemState.bodies.length);
    }
    let bodiesToRender: Body[] = systemState.bodies;
    let mode = this.state.cameraMode;
    if (mode === CameraMode.SHIP_RELATIVE) {
      const playerShipId = this.state.playerShipId;
      const playerShip = playerShipId ? systemState.bodies.find(b => b.id === playerShipId) as Ship | undefined : undefined;
      this.cameraManager.update(mode, { playerShip });
      bodiesToRender = systemState.bodies.filter(b => b.id !== this.state.playerShipId);
    } else if (mode !== CameraMode.SYSTEM_ORBITAL) {
      // Galaxy camera: build view-projection and extract right/up for billboards
      const proj = mat4.create();
      const aspect = this.textureSize.height > 0 ? (this.textureSize.width / this.textureSize.height) : 16/9;
      mat4.perspective(proj, Math.PI / 4, aspect, 0.1, 2000.0);
      const view = mat4.create();
      // Use camera.look_at as target (class uses look_at)
      mat4.lookAt(view, this.camera.eye as unknown as number[], this.camera.look_at as unknown as number[], this.camera.up as unknown as number[]);
      const viewProj = mat4.multiply(mat4.create(), proj, view);

      // Extract camera right/up from view matrix columns
      const camRight: Vec3 = [view[0], view[4], view[8]];
      const camUp: Vec3 = [view[1], view[5], view[9]];

      const camData = new Float32Array(32); // 128 bytes
      camData.set(viewProj, 0);
      camData.set(camRight, 16);
      camData.set(camUp, 20);
      this.device.queue.writeBuffer(this.galaxyCameraUniformBuffer, 0, camData);
    this.setTheme(this.currentThemeName, this.currentResponseName, deltaTime);
    this._render();
    // HUD draw for galactic map: currently no bodies data; skip
    this.input.tick();
      requestAnimationFrame(() => this._update());
      return;
    }

    // Determine render scale: fit system in SYSTEM_ORBITAL, otherwise 1:1
    let currentRenderScale: number;
    if (mode === CameraMode.SYSTEM_ORBITAL) {
        const systemRadius = systemState.bodies.reduce((maxDist, b) => {
            const dist = Math.hypot(b.position[0], b.position[1], b.position[2]);
            return Math.max(maxDist, dist);
        }, 0);
        currentRenderScale = systemRadius > 0 ? VISUAL_SETTINGS.systemViewSize / systemRadius : 1.0;
    } else {
        currentRenderScale = 1.0;
    }

    // Update camera by controller for SYSTEM_ORBITAL prior to serialization
    if (mode === CameraMode.SYSTEM_ORBITAL) {
      this.cameraManager.update(mode, {
        bodies: systemState.bodies,
        scale: currentRenderScale,
        viewport: { width: this.textureSize.width, height: this.textureSize.height },
        vfov: 25.0,
      });
    }
    // Serialize with selected render scale
    const tmp = this._serializeSystemState(bodiesToRender, this.camera.focusBodyId, currentRenderScale);
    let sphereData: Float32Array = tmp.sphereData;
    let systemScale: number = tmp.systemScale;

    if (mode === CameraMode.SYSTEM_ORBITAL) {
      this.lastSystemScale = systemScale;
    } else {
      this.lastSystemScale = systemScale;
    }
    // Upload camera uniforms. The world is camera-relative, so the camera itself is at the origin.
    const cameraData = new Float32Array(16);
    const cam = this.camera;
    cameraData.set([0, 0, 0], 0);
    const lookAtRelative: Vec3 = [
      cam.look_at[0] - cam.eye[0],
      cam.look_at[1] - cam.eye[1],
      cam.look_at[2] - cam.eye[2],
    ];
    cameraData.set(lookAtRelative, 4);
    cameraData.set(cam.up, 8);
    cameraData[12] = 0.0;
    this.device.queue.writeBuffer(this.cameraUniformBuffer, 0, cameraData);

    // Pad sphere buffer to full body count to match compute shader expectations
    const floatsPerSphere = 16;
    const expectedCount = this.lastKnownBodyCount;
    const upload = new Float32Array(Math.max(1, expectedCount * floatsPerSphere));
    upload.set(sphereData, 0);
    this.device.queue.writeBuffer(this.spheresBuffer, 0, upload.buffer as ArrayBuffer, 0, upload.byteLength);
    this._lastBodiesForMask = systemState.bodies;

    // Update analytic orbits and uniforms if enabled
    if (this.state.showOrbits) {
      this._updateAnalyticOrbits(systemState.bodies, systemScale);
      const proj = mat4.create();
      const aspect = this.textureSize.height > 0 ? (this.textureSize.width / this.textureSize.height) : 16/9;
      mat4.perspective(proj, Math.PI * 25 / 180, aspect, 0.00001, 10000.0);
      const view = mat4.create();
      mat4.lookAt(view, this.camera.eye as unknown as number[], this.camera.look_at as unknown as number[], this.camera.up as unknown as number[]);
      const viewProj = mat4.multiply(mat4.create(), proj, view);
      const theme = themes[this.currentThemeName as keyof typeof themes];
      const color = theme.accent;
      const u = new Float32Array(20);
      u.set(viewProj, 0);
      u.set(color, 16);
      this.device.queue.writeBuffer(this.orbitsUniformBuffer, 0, u);
    }

    this.setTheme(this.currentThemeName, this.currentResponseName, deltaTime);
    this._render();
    // Draw HUD markers
    if (this.hud) {
      const viewport = { width: this.textureSize.width, height: this.textureSize.height };
      if (mode === CameraMode.SYSTEM_ORBITAL) {
        const scaledBodies = systemState.bodies.map(b => ({
          ...b,
          position: [
            b.position[0] * currentRenderScale,
            b.position[1] * currentRenderScale,
            b.position[2] * currentRenderScale,
          ] as Vec3,
        }));
        this.hud.draw(scaledBodies, this.camera, viewport, mode, this.state.playerShipId);
      } else {
        this.hud.draw(systemState.bodies, this.camera, viewport, mode, this.state.playerShipId);
      }
    }
    this.input.tick();
    requestAnimationFrame(() => this._update());
  }

  private _render() {
    const commandEncoder = this.device.createCommandEncoder();
    const themeBg = themes[this.currentThemeName as keyof typeof themes]?.bg ?? [0, 0, 0];

    if (this.state.cameraMode !== CameraMode.GALACTIC_MAP) {
      const sourceTexture = this.frameCount % 2 === 0 ? this.postFxTextureB : this.postFxTextureA;
      const destinationTexture = this.frameCount % 2 === 0 ? this.postFxTextureA : this.postFxTextureB;

      // Build target info uniform: only the focused body position as a single target
      {
        const buf = new ArrayBuffer(272);
        const u32 = new Uint32Array(buf);
        const f32 = new Float32Array(buf);
        u32[0] = 1; // count = 1
        const radiusPx = 20.0;
        f32[1] = radiusPx;
        f32[2] = this.textureSize.width;
        f32[3] = this.textureSize.height;
        const scale = this.lastSystemScale;
        const focusId = this.camera.focusBodyId;
        const b = this._lastBodiesForMask.find(bb => bb.id === focusId) || this._lastBodiesForMask[0];
        if (b) {
          const x = (b.position[0] * scale) - this.camera.eye[0];
          const y = (b.position[1] * scale) - this.camera.eye[1];
          const z = (b.position[2] * scale) - this.camera.eye[2];
          const look = [this.camera.look_at[0] - this.camera.eye[0], this.camera.look_at[1] - this.camera.eye[1], this.camera.look_at[2] - this.camera.eye[2]] as unknown as number[];
          const forwardLen = Math.hypot(look[0], look[1], look[2]) || 1;
          const fx = look[0] / forwardLen, fy = look[1] / forwardLen, fz = look[2] / forwardLen;
          const up = this.camera.up as unknown as number[];
          const rx = up[1]*fz - up[2]*fy;
          const ry = up[2]*fx - up[0]*fz;
          const rz = up[0]*fy - up[1]*fx;
          const rLen = Math.hypot(rx, ry, rz) || 1;
          const rnx = rx/rLen, rny = ry/rLen, rnz = rz/rLen;
          const ux = fy*rnz - fz*rny;
          const uy = fz*rnx - fx*rnz;
          const uz = fx*rny - fy*rnx;
          const vx = x*rnx + y*rny + z*rnz;
          const vy = x*ux + y*uy + z*uz;
          const vz = x*fx + y*fy + z*fz;
          const fov = 25 * Math.PI / 180;
          const projScale = 1 / Math.tan(fov/2);
          const ndcX = (vx / Math.max(1e-4, vz)) * projScale;
          const ndcY = (vy / Math.max(1e-4, vz)) * projScale;
          const px = (ndcX * 0.5 + 0.5) * this.textureSize.width;
          const py = (1.0 - (ndcY * 0.5 + 0.5)) * this.textureSize.height;
          const base = 4;
          f32[base+0] = px;
          f32[base+1] = py;
          f32[base+2] = 0.0;
          f32[base+3] = 0.0;
        }
        this.device.queue.writeBuffer(this.targetInfoUniform, 0, buf);
      }

      const postFxBindGroup = this.device.createBindGroup({
        layout: this.renderPipeline.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: this.sampler },
          { binding: 1, resource: this.computeShader.texture.createView() },
          { binding: 2, resource: { buffer: this.themeUniformBuffer } },
          { binding: 3, resource: sourceTexture.createView() },
          { binding: 4, resource: this.orbitsTexture.createView() },
          { binding: 5, resource: { buffer: this.targetInfoUniform } },
        ]
      });

      const computePass = commandEncoder.beginComputePass();
      computePass.setPipeline(this.computeShader.pipeline);
      computePass.setBindGroup(0, this.computeShader.bindGroup);
      computePass.dispatchWorkgroups(Math.ceil(this.textureSize.width / 8), Math.ceil(this.textureSize.height / 8));
      computePass.end();

      // Pass 1: Render PostFX into offscreen destination ping-pong texture
      const postPass = commandEncoder.beginRenderPass({
        colorAttachments: [{
          view: destinationTexture.createView(),
          loadOp: "clear",
          storeOp: "store",
          clearValue: { r: themeBg[0], g: themeBg[1], b: themeBg[2], a: 1 }
        }]
      });
      postPass.setPipeline(this.renderPipeline);
      postPass.setBindGroup(0, postFxBindGroup);
      postPass.draw(6);
      postPass.end();

      // Optional Pass 1.5: Draw orbits into the orbits texture, then composite in PostFX
      if (this.state.showOrbits) {
        const orbitsPass = commandEncoder.beginRenderPass({
          colorAttachments: [{
            view: this.orbitsTexture.createView(),
            loadOp: 'clear',
            storeOp: 'store',
            clearValue: { r: 0, g: 0, b: 0, a: 1 },
          }]
        });
        orbitsPass.setPipeline(this.orbitsPipeline);
        orbitsPass.setBindGroup(0, this.orbitsBindGroup);
        for (let i = 0; i < this.orbitBuffers.length; i++) {
          const count = Math.max(0, this.orbitCounts[i] - this.ORBIT_TAIL_HIDE);
          if (count < 2) continue;
          orbitsPass.setVertexBuffer(0, this.orbitBuffers[i]);
          orbitsPass.draw(count);
        }
        orbitsPass.end();
      }

      // Pass 2: Present the destination texture to the canvas
      const currentTexture = this.context.getCurrentTexture();
      const presentPass = commandEncoder.beginRenderPass({
        colorAttachments: [{
          view: currentTexture.createView(),
          loadOp: "clear",
          storeOp: "store",
          clearValue: { r: themeBg[0], g: themeBg[1], b: themeBg[2], a: 1 }
        }]
      });
      // Build a bind group that feeds the destination texture into the present shader at binding 1
      const presentBindGroup = this.device.createBindGroup({
        layout: this.presentPipeline.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: this.sampler },
          { binding: 1, resource: destinationTexture.createView() },
        ]
      });
      presentPass.setPipeline(this.presentPipeline);
      presentPass.setBindGroup(0, presentBindGroup);
      presentPass.draw(6);
      presentPass.end();

    } else {
      const currentTexture = this.context.getCurrentTexture();
      const renderPass = commandEncoder.beginRenderPass({
        colorAttachments: [{
          view: currentTexture.createView(),
          loadOp: "clear",
          storeOp: "store",
          clearValue: { r: 0.08, g: 0.0, b: 0.12, a: 1.0 } // purple background for galaxy view
        }],
        depthStencilAttachment: {
          view: this.galaxyDepthTexture.createView(),
          depthClearValue: 1.0,
          depthLoadOp: 'clear',
          depthStoreOp: 'store',
        },
      });
      renderPass.setPipeline(this.galaxyPipeline);
      renderPass.setBindGroup(0, this.galaxyBindGroup);
      renderPass.draw(4, this.galaxy.stars.length);
      renderPass.end();
    }

    this.device.queue.submit([commandEncoder.finish()]);
    this.frameCount++;
  }

  private _createGalaxyPipeline() {
    const module = this.device.createShaderModule({ code: galaxyShaderWGSL });
    this.galaxyPipeline = this.device.createRenderPipeline({
      layout: "auto",
      vertex: { module, entryPoint: "vertexMain" },
      fragment: { module, entryPoint: "fragmentMain", targets: [{ format: this.presentationFormat }] },
      primitive: { topology: "triangle-strip" },
      depthStencil: {
        format: 'depth24plus',
        depthWriteEnabled: true,
        depthCompare: 'less',
      },
    });
  }

  private _createOrbitsPipeline() {
    const module = this.device.createShaderModule({ code: orbitsShaderWGSL });
    this.orbitsPipeline = this.device.createRenderPipeline({
      layout: 'auto',
      vertex: {
        module,
        entryPoint: 'vertexMain',
        buffers: [
          {
            arrayStride: 12,
            attributes: [
              { shaderLocation: 0, format: 'float32x3', offset: 0 },
            ],
          },
        ],
      },
      fragment: {
        module,
        entryPoint: 'fragmentMain',
        targets: [
          {
            format: this.presentationFormat,
            // Additive blend to overlay lines
            blend: {
              color: { srcFactor: 'one', dstFactor: 'one', operation: 'add' },
              alpha: { srcFactor: 'one', dstFactor: 'one', operation: 'add' },
            },
          },
        ],
      },
      primitive: { topology: 'line-strip' },
    });
  }

  // Compute analytic Keplerian ellipse (osculating) for each body around its strongest parent
  private _updateAnalyticOrbits(bodies: Body[], systemScale: number) {
    const hierarchy = buildSystemHierarchy(bodies);
    for (let i = 0; i < bodies.length; i++) {
      const body = bodies[i];
      const parentId = hierarchy.get(body.id);
      if (!parentId) {
        this.orbitCounts[i] = 0;
        continue;
      }
      const parent = bodies.find(b => b.id === parentId);
      if (!parent) {
        this.orbitCounts[i] = 0;
        continue;
      }

      // Relative state in world km
      const rx = body.position[0] - parent.position[0];
      const ry = body.position[1] - parent.position[1];
      const rz = body.position[2] - parent.position[2];
      const vx = body.velocity[0] - parent.velocity[0];
      const vy = body.velocity[1] - parent.velocity[1];
      const vz = body.velocity[2] - parent.velocity[2];

      const r = Math.hypot(rx, ry, rz);
      const v = Math.hypot(vx, vy, vz);
      const mu = G * (parent.mass + body.mass);

      // Specific angular momentum h = r x v
      const hx = ry * vz - rz * vy;
      const hy = rz * vx - rx * vz;
      const hz = rx * vy - ry * vx;
      const h = Math.hypot(hx, hy, hz) || 1e-9;

      // Eccentricity vector e = (v x h)/mu - r/|r|
      const cx = vy * hz - vz * hy;
      const cy = vz * hx - vx * hz;
      const cz = vx * hy - vy * hx;
      const ex = cx / mu - rx / r;
      const ey = cy / mu - ry / r;
      const ez = cz / mu - rz / r;
      const e = Math.hypot(ex, ey, ez);

      // Semi-major axis a from vis-viva: v^2 = mu*(2/r - 1/a)
      const invA = 2 / r - (v * v) / mu;
      if (!isFinite(invA) || invA <= 0) {
        this.orbitCounts[i] = 0; // Parabolic/hyperbolic: skip
        continue;
      }
      const a = 1 / invA;
      const bSemi = a * Math.sqrt(Math.max(0, 1 - e * e));

      // Perifocal basis (P, Q, W)
      const hnx = hx / h, hny = hy / h, hnz = hz / h; // W = h-hat
      // P along eccentricity (periapsis direction). If nearly circular, use radial direction.
      let Px: number, Py: number, Pz: number;
      if (e > 1e-6) {
        const invEL = 1.0 / e;
        Px = ex * invEL; Py = ey * invEL; Pz = ez * invEL;
      } else {
        const invR = 1.0 / Math.max(1e-9, r);
        Px = rx * invR; Py = ry * invR; Pz = rz * invR;
      }
      // Ensure P is orthogonalized within the plane (remove any tiny component along W)
      const dotPW = Px * hnx + Py * hny + Pz * hnz;
      Px -= dotPW * hnx; Py -= dotPW * hny; Pz -= dotPW * hnz;
      const pLen = Math.hypot(Px, Py, Pz) || 1;
      Px /= pLen; Py /= pLen; Pz /= pLen;
      // Q = W x P
      let Qx = hny * Pz - hnz * Py;
      let Qy = hnz * Px - hnx * Pz;
      let Qz = hnx * Py - hny * Px;
      const qLen = Math.hypot(Qx, Qy, Qz) || 1;
      Qx /= qLen; Qy /= qLen; Qz /= qLen;

      // Center at parent position; generate ellipse points in world km, then scale to render units
      const arr = this.orbitCPUData[i];
      const N = Math.min(this.ORBIT_SAMPLES + 1, this.ORBIT_MAX_POINTS); // close the ring by repeating first point
      for (let k = 0; k < N; k++) {
        const theta = (k / (N - 1)) * Math.PI * 2;
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        const xPeri = a * (cosT - e);
        const yPeri = bSemi * sinT;
        const wx = Px * xPeri + Qx * yPeri;
        const wy = Py * xPeri + Qy * yPeri;
        const wz = Pz * xPeri + Qz * yPeri;
        const base = k * 3;
        arr[base + 0] = (parent.position[0] + wx) * systemScale;
        arr[base + 1] = (parent.position[1] + wy) * systemScale;
        arr[base + 2] = (parent.position[2] + wz) * systemScale;
      }
      this.orbitCounts[i] = N;
      const bytes = N * 3 * 4;
      this.device.queue.writeBuffer(this.orbitBuffers[i], 0, arr.buffer as ArrayBuffer, 0, bytes);
    }
  }

  public clearOrbitHistory(): void {
    for (let i = 0; i < this.orbitCounts.length; i++) {
      this.orbitCounts[i] = 0;
    }
  }

  private _serializeStars(stars: Star[]): Float32Array {
    // Layout: position.xyz, color.xyz, size (with implicit alignment)
    const floatsPerStar = 8;
    const data = new Float32Array(stars.length * floatsPerStar);
    for (let i = 0; i < stars.length; i++) {
      const base = i * floatsPerStar;
      const s = stars[i];
      data.set(s.position, base + 0);
      data.set(s.color, base + 4);
      data[base + 7] = s.size;
    }
    return data;
  }
}
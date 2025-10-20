import computeShaderWGSL from './shaders/compute.wgsl?raw';
import postfxShaderWGSL from './shaders/postfx.wgsl?raw';
import galaxyShaderWGSL from './shaders/galaxy.wgsl?raw';
import { AppState, ViewMode } from './state';
import orbitsShaderWGSL from './shaders/orbits.wgsl?raw';
import { themes } from './theme';
import { spectralResponses } from './spectral';
import type { Authority } from './authority/authority';
import type { Body } from './shared/types';
import { G } from './shared/constants';
import { Camera } from './camera';
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

export class Renderer {

  private canvas: HTMLCanvasElement;
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
  public camera!: Camera;
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
  // Orbits overlay
  private orbitsPipeline!: GPURenderPipeline;
  private orbitsUniformBuffer!: GPUBuffer;
  private orbitsBindGroup!: GPUBindGroup;
  private orbitBuffers: GPUBuffer[] = [];
  private orbitCPUData: Float32Array[] = [];
  private orbitCounts: number[] = [];
  private readonly ORBIT_MAX_POINTS = 2048;

  constructor(canvas: HTMLCanvasElement, authority: Authority, state: AppState) {
    this.canvas = canvas;
    this.authority = authority;
    this.camera = new Camera(this.canvas);
    this.state = state;
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

    // Initialize camera uniforms before any rendering
    {
      const cameraData = new Float32Array(16);
      const eye = this.camera.eye;
      const lookAt = this.camera.look_at;
      const dist = Math.hypot(eye[0] - lookAt[0], eye[1] - lookAt[1], eye[2] - lookAt[2]);
      // Camera-relative uniforms: eye at origin, look_at relative to eye
      cameraData.set([0, 0, 0], 0);
      const relativeLookAt: Vec3 = [
        lookAt[0] - eye[0],
        lookAt[1] - eye[1],
        lookAt[2] - eye[2],
      ];
      cameraData.set(relativeLookAt, 4);
      cameraData.set(this.camera.up, 8);
      cameraData[12] = dist; // distance_to_target
      this.device.queue.writeBuffer(this.cameraUniformBuffer, 0, cameraData);
    }

    const width = this.canvas.getBoundingClientRect().width;
    const aspectRation = 16.0 / 9.0;
    const height = Math.floor(width / aspectRation);
    this.canvas.width = Math.max(1, Math.min(width, this.device.limits.maxTextureDimension2D));
    this.canvas.height = Math.max(1, Math.min(height, this.device.limits.maxTextureDimension2D));
    this.textureSize = { width: this.canvas.width, height: this.canvas.height };

    const systemState = await this.authority.query();
    const initSerialized = this._serializeSystemState(systemState.bodies, this.camera.focusBodyId);
    const sphereData = initSerialized.sphereData;
    const numBodies = systemState.bodies.length;

    this.spheresBuffer = this.device.createBuffer({
      size: sphereData.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(this.spheresBuffer, 0, sphereData.buffer as ArrayBuffer, 0, sphereData.byteLength);

    this.computeShader = this.createComputeShader(this.textureSize, numBodies, this.spheresBuffer);
    this.renderPipeline = this.createPostFXPipeline(this.presentationFormat);
    this.presentPipeline = this.createPresentPipeline(this.presentationFormat);

    // Orbits pipeline & buffers
    this._createOrbitsPipeline();
    // 80 bytes for alignment (64 for viewProjection + 16 padding, store color in first 3 of next vec4)
    this.orbitsUniformBuffer = this.device.createBuffer({
      size: 80,
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
        // Recreate galaxy depth texture on resize
        this.galaxyDepthTexture = this.device.createTexture({
          size: this.textureSize,
          format: 'depth24plus',
          usage: GPUTextureUsage.RENDER_ATTACHMENT,
        });
        this._render();
        console.log("Rerender time:", performance.now() - renderTime);
      }
    });
    observer.observe(this.canvas);

    this._update();
  }

  private _serializeSystemState(bodies: Body[], focusBodyId: string | null): { sphereData: Float32Array, focusBodyRenderedRadius: number } {
    const floatsPerSphere = 16; // 64 bytes per sphere to match WGSL layout
    const sphereData = new Float32Array(bodies.length * floatsPerSphere);
    const sphereDataU32 = new Uint32Array(sphereData.buffer);
    const scale = 1 / 1.5e8; // Scales 1 AU to 1 unit
    const eye: Vec3 = this.camera.eye;

    const SYSTEM_EXAGGERATION = 100.0;
    let focusBodyRenderedRadius = 0;

    const focusBody = bodies.find(b => b.id === focusBodyId) || null;
    const centralBody = bodies.reduce((a, b) => a.mass > b.mass ? a : b);
    const contextParent = focusBody && focusBody.id !== centralBody.id ? focusBody : centralBody;

    // Determine the exaggeration factor for the current context
    let contextExaggeration = SYSTEM_EXAGGERATION;
    if (contextParent !== centralBody) {
      // Consider bodies less massive than the context parent and not the central body as potential moons
      const moons = bodies.filter(b => b.mass < contextParent.mass && b.id !== centralBody.id);
      if (moons.length > 0) {
        const closestMoonDist = Math.min(...moons.map(moon => {
          const dx = moon.position[0] - contextParent.position[0];
          const dy = moon.position[1] - contextParent.position[1];
          const dz = moon.position[2] - contextParent.position[2];
          return Math.sqrt(dx*dx + dy*dy + dz*dz);
        }));
        const safetyMargin = 0.8;
        const maxSafeRadius = (closestMoonDist * scale) * safetyMargin;
        const baseRadius = contextParent.radius * scale;
        contextExaggeration = maxSafeRadius / baseRadius;
      }
    }

    bodies.forEach((body, i) => {
      const f_base = i * floatsPerSphere; // float index base
      const u_base = f_base;              // u32 index base

      // Exaggeration factor per-body based on context
      let exaggerationFactor = SYSTEM_EXAGGERATION;
      if (contextParent !== centralBody) {
        const isContext = body.id === contextParent.id;
        const isPotentialMoon = body.mass < contextParent.mass && body.id !== centralBody.id;
        if (isContext || isPotentialMoon) {
          exaggerationFactor = contextExaggeration;
        }
      }

      const renderedRadius = body.radius * scale * exaggerationFactor;
      if (body.id === focusBodyId) {
        focusBodyRenderedRadius = renderedRadius;
      }

      // Sphere.center (vec3f) in camera-relative coordinates -> slots 0, 1, 2
      sphereData[f_base + 0] = (body.position[0] * scale) - eye[0];
      sphereData[f_base + 1] = (body.position[1] * scale) - eye[1];
      sphereData[f_base + 2] = (body.position[2] * scale) - eye[2];

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
    return { sphereData, focusBodyRenderedRadius };
  }

  private async _update() {
    const deltaTime = 1 / 60.0;
    await this.authority.tick(deltaTime);
    const systemState = await this.authority.query();
    if (this.state.viewMode === ViewMode.System) {
      // Update camera target based on focused body (scale matches serialize)
      this.camera.update(systemState.bodies, 1 / 1.5e8);
      const { sphereData, focusBodyRenderedRadius } = this._serializeSystemState(systemState.bodies, this.camera.focusBodyId);
      if (this.camera.pendingFrame) {
        this.camera.completePendingFrame(focusBodyRenderedRadius);
      }
      // Upload camera uniforms (eye, target, up, distance_to_target)
      const cameraData = new Float32Array(16);
      const eye = this.camera.eye;
      const lookAt = this.camera.look_at;
      const dist = Math.hypot(eye[0] - lookAt[0], eye[1] - lookAt[1], eye[2] - lookAt[2]);
      // Camera-relative uniforms: eye at origin, look_at relative to eye
      cameraData.set([0, 0, 0], 0);
      const relativeLookAt: Vec3 = [
        lookAt[0] - eye[0],
        lookAt[1] - eye[1],
        lookAt[2] - eye[2],
      ];
      cameraData.set(relativeLookAt, 4);
      cameraData.set(this.camera.up, 8);
      cameraData[12] = dist; // distance_to_target
      this.device.queue.writeBuffer(this.cameraUniformBuffer, 0, cameraData);
      this.device.queue.writeBuffer(this.spheresBuffer, 0, sphereData.buffer as ArrayBuffer, 0, sphereData.byteLength);

      // Update orbit trails and uniforms if enabled
      if (this.state.showOrbits) {
        const scale = 1 / 1.5e8;
        // Append current positions (scaled) to trails
        for (let i = 0; i < systemState.bodies.length; i++) {
          const b = systemState.bodies[i];
          const arr = this.orbitCPUData[i];
          let n = this.orbitCounts[i];
          const px = b.position[0] * scale;
          const py = b.position[1] * scale;
          const pz = b.position[2] * scale;
          if (n < this.ORBIT_MAX_POINTS) {
            const base = n * 3;
            arr[base + 0] = px;
            arr[base + 1] = py;
            arr[base + 2] = pz;
            n++;
          } else {
            // shift left by one vec3
            arr.copyWithin(0, 3);
            const base = (this.ORBIT_MAX_POINTS - 1) * 3;
            arr[base + 0] = px;
            arr[base + 1] = py;
            arr[base + 2] = pz;
            n = this.ORBIT_MAX_POINTS;
          }
          this.orbitCounts[i] = n;
          // upload current trail
          const bytes = n * 3 * 4;
          this.device.queue.writeBuffer(this.orbitBuffers[i], 0, arr.buffer as ArrayBuffer, 0, bytes);
        }

        // Update orbits uniforms (viewProjection + color)
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
        u.set(color, 16); // write into vec4 slot starting at 16 (only first 3 used)
        this.device.queue.writeBuffer(this.orbitsUniformBuffer, 0, u);
      }
    } else {
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
    }
    this.setTheme(this.currentThemeName, this.currentResponseName, deltaTime);
    this._render();
    requestAnimationFrame(() => this._update());
  }

  private _render() {
    const commandEncoder = this.device.createCommandEncoder();

    if (this.state.viewMode === ViewMode.System) {
      const sourceTexture = this.frameCount % 2 === 0 ? this.postFxTextureB : this.postFxTextureA;
      const destinationTexture = this.frameCount % 2 === 0 ? this.postFxTextureA : this.postFxTextureB;

      // Create per-frame bind group to include prev frame texture at binding 3
      const postFxBindGroup = this.device.createBindGroup({
        layout: this.renderPipeline.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: this.sampler },
          { binding: 1, resource: this.computeShader.texture.createView() },
          { binding: 2, resource: { buffer: this.themeUniformBuffer } },
          { binding: 3, resource: sourceTexture.createView() },
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
          clearValue: { r: 0, g: 0, b: 0, a: 1 }
        }]
      });
      postPass.setPipeline(this.renderPipeline);
      postPass.setBindGroup(0, postFxBindGroup);
      postPass.draw(6);
      postPass.end();

      // Optional Pass 1.5: Draw orbits into the destination texture
      if (this.state.showOrbits) {
        const orbitsPass = commandEncoder.beginRenderPass({
          colorAttachments: [{
            view: destinationTexture.createView(),
            loadOp: 'load',
            storeOp: 'store',
          }]
        });
        orbitsPass.setPipeline(this.orbitsPipeline);
        orbitsPass.setBindGroup(0, this.orbitsBindGroup);
        for (let i = 0; i < this.orbitBuffers.length; i++) {
          const count = this.orbitCounts[i];
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
          clearValue: { r: 0, g: 0, b: 0, a: 1 }
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
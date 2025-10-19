import computeShaderWGSL from './shaders/compute.wgsl?raw';
import postfxShaderWGSL from './shaders/postfx.wgsl?raw';
import galaxyShaderWGSL from './shaders/galaxy.wgsl?raw';
import { AppState, ViewMode } from './state';
import { themes } from './theme';
import { spectralResponses } from './spectral';
import type { Authority } from './authority/authority';
import type { Body } from './shared/types';
import { Camera } from './camera';
import { Galaxy } from './galaxy';
import type { Star } from './galaxy';
import type { Vec3 } from './shared/types';
import { mat4 } from 'gl-matrix';

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
      const cameraData = new Float32Array(12);
      cameraData.set(this.camera.eye, 0);
      cameraData.set(this.camera.look_at, 4);
      cameraData.set(this.camera.up, 8);
      this.device.queue.writeBuffer(this.cameraUniformBuffer, 0, cameraData);
    }

    const width = this.canvas.getBoundingClientRect().width;
    const aspectRation = 16.0 / 9.0;
    const height = Math.floor(width / aspectRation);
    this.canvas.width = Math.max(1, Math.min(width, this.device.limits.maxTextureDimension2D));
    this.canvas.height = Math.max(1, Math.min(height, this.device.limits.maxTextureDimension2D));
    this.textureSize = { width: this.canvas.width, height: this.canvas.height };

    const systemState = await this.authority.query();
    const sphereData = this._serializeSystemState(systemState.bodies);
    const numBodies = systemState.bodies.length;

    this.spheresBuffer = this.device.createBuffer({
      size: sphereData.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(this.spheresBuffer, 0, sphereData.buffer as ArrayBuffer, 0, sphereData.byteLength);

    this.computeShader = this.createComputeShader(this.textureSize, numBodies, this.spheresBuffer);
    this.renderPipeline = this.createPostFXPipeline(this.presentationFormat);
    this.presentPipeline = this.createPresentPipeline(this.presentationFormat);

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

    const observer = new ResizeObserver(entries => {
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
        this.computeShader = this.createComputeShader(this.textureSize, numBodies, this.spheresBuffer);
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

  private _serializeSystemState(bodies: Body[]): Float32Array {
    const floatsPerSphere = 16; // 64 bytes per sphere to match WGSL layout
    const sphereData = new Float32Array(bodies.length * floatsPerSphere);
    const sphereDataU32 = new Uint32Array(sphereData.buffer);

    bodies.forEach((body, i) => {
      const f_base = i * floatsPerSphere; // float index base
      const u_base = f_base;              // u32 index base
      const scale = 1 / 1.5e8; // Scales 1 AU to 1 unit

      // Sphere.center (vec3f) -> slots 0, 1, 2
      sphereData[f_base + 0] = body.position[0] * scale;
      sphereData[f_base + 1] = body.position[1] * scale;
      sphereData[f_base + 2] = body.position[2] * scale - 2.0;

      // Sphere.radius (f32) -> slot 3
      sphereData[f_base + 3] = body.radius * scale * 100;

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
    return sphereData;
  }

  private async _update() {
    const deltaTime = 1 / 60.0;
    await this.authority.tick(deltaTime);
    const systemState = await this.authority.query();
    if (this.state.viewMode === ViewMode.System) {
      // Update camera target based on focused body (scale matches serialize)
      this.camera.update(systemState.bodies, 1 / 1.5e8);
      const sphereData = this._serializeSystemState(systemState.bodies);
      // Upload camera uniforms (eye, target, up)
      const cameraData = new Float32Array(12);
      cameraData.set(this.camera.eye, 0);
      cameraData.set(this.camera.look_at, 4);
      cameraData.set(this.camera.up, 8);
      this.device.queue.writeBuffer(this.cameraUniformBuffer, 0, cameraData);
      this.device.queue.writeBuffer(this.spheresBuffer, 0, sphereData.buffer as ArrayBuffer, 0, sphereData.byteLength);
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
    this.setTheme('amber', 'Visible (Y)', deltaTime);
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
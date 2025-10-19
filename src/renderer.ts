import computeShaderWGSL from './shaders/compute.wgsl?raw';
import renderShaderWGSL from './shaders/render.wgsl?raw';
import type { Authority } from './authority/authority';
import type { Body } from './shared/types';
import { Camera } from './camera';

export class Renderer {

  private canvas: HTMLCanvasElement;
  private device!: GPUDevice;
  private context!: GPUCanvasContext;
  private presentationFormat!: GPUTextureFormat;
  private textureSize!: { width: number, height: number };
  private computeShader!: { pipeline: GPUComputePipeline, bindGroup: GPUBindGroup, texture: GPUTexture };
  private renderPipeline!: GPURenderPipeline;
  private sampler!: GPUSampler;
  private renderBindGroup!: GPUBindGroup;
  private authority: Authority;
  private spheresBuffer!: GPUBuffer;
  private cameraUniformBuffer!: GPUBuffer;
  public camera!: Camera;

  constructor(canvas: HTMLCanvasElement, authority: Authority) {
    this.canvas = canvas;
    this.authority = authority;
    this.camera = new Camera(this.canvas);
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
      format: 'rgba8unorm',
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

  private createRenderPipeline(format: GPUTextureFormat) {
    const module = this.device.createShaderModule({
      label: "Render shader",
      code: renderShaderWGSL
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
        topology: "triangle-strip",
        stripIndexFormat: "uint32"
      }
    });
  }

  private updateBindGroups() {
    this.renderBindGroup = this.device.createBindGroup({
      layout: this.renderPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: this.sampler },
        { binding: 1, resource: this.computeShader.texture.createView() }
      ]
    });
  }

  public async start() {
    const { device, context, presentationFormat } = await this.initWebGPU();
    this.device = device;
    this.context = context;
    this.presentationFormat = presentationFormat;

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
    this.renderPipeline = this.createRenderPipeline(this.presentationFormat);

    this.sampler = this.device.createSampler({
      magFilter: "linear",
      minFilter: "linear"
    });

    this.updateBindGroups();

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
        this.updateBindGroups();
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
    await this.authority.tick(1 / 60.0);
    const systemState = await this.authority.query();
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
    this._render();
    requestAnimationFrame(() => this._update());
  }

  private _render() {
    const commandEncoder = this.device.createCommandEncoder();

    const computePass = commandEncoder.beginComputePass();
    computePass.setPipeline(this.computeShader.pipeline);
    computePass.setBindGroup(0, this.computeShader.bindGroup);
    computePass.dispatchWorkgroups(Math.ceil(this.textureSize.width / 8), Math.ceil(this.textureSize.height / 8));
    computePass.end();

    const renderPass = commandEncoder.beginRenderPass({
      colorAttachments: [{
        view: this.context.getCurrentTexture().createView(),
        loadOp: "clear",
        storeOp: "store",
        clearValue: { r: 0, g: 0, b: 0, a: 1 }
      }]
    });
    renderPass.setPipeline(this.renderPipeline);
    renderPass.setBindGroup(0, this.renderBindGroup);
    renderPass.draw(4);
    renderPass.end();

    this.device.queue.submit([commandEncoder.finish()]);
  }
}
import computeShaderWGSL from './shaders/compute.wgsl?raw';
import renderShaderWGSL from './shaders/render.wgsl?raw';
import { createRandomScene, serializeScene } from './scene';

export class Renderer {
  private static readonly NUM_SPHERES = 30;

  private canvas: HTMLCanvasElement;
  private device!: GPUDevice;
  private context!: GPUCanvasContext;
  private presentationFormat!: GPUTextureFormat;
  private textureSize!: { width: number, height: number };
  private computeShader!: { pipeline: GPUComputePipeline, bindGroup: GPUBindGroup, texture: GPUTexture };
  private renderPipeline!: GPURenderPipeline;
  private sampler!: GPUSampler;
  private renderBindGroup!: GPUBindGroup;

  constructor(canvas: HTMLCanvasElement) {
    this.canvas = canvas;
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

  private createComputeShader(textureSize: { width: number, height: number }) {
    const module = this.device.createShaderModule({
      label: "Compute shader",
      code: computeShaderWGSL
    });

    const pipeline = this.device.createComputePipeline({
      label: "Compute pipeline",
      layout: "auto",
      compute: {
        module,
        entryPoint: "main"
      }
    });

    const spheres = createRandomScene();
    const sphereData = serializeScene(spheres);

    const spheresBuffer = this.device.createBuffer({
      size: sphereData.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

    // Use ArrayBuffer overload to satisfy types
    this.device.queue.writeBuffer(spheresBuffer, 0, sphereData.buffer as ArrayBuffer, 0, sphereData.byteLength);

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
        { binding: 1, resource: texture.createView() }
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

  private render() {
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

  public async start() {
    const { device, context, presentationFormat } = await this.initWebGPU();
    this.device = device;
    this.context = context;
    this.presentationFormat = presentationFormat;

    const width = this.canvas.getBoundingClientRect().width;
    const aspectRation = 16.0 / 9.0;
    const height = Math.floor(width / aspectRation);
    this.canvas.width = Math.max(1, Math.min(width, this.device.limits.maxTextureDimension2D));
    this.canvas.height = Math.max(1, Math.min(height, this.device.limits.maxTextureDimension2D));
    this.textureSize = { width: this.canvas.width, height: this.canvas.height };

    this.computeShader = this.createComputeShader(this.textureSize);
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
        this.computeShader = this.createComputeShader(this.textureSize);
        this.updateBindGroups();
        this.render();
        console.log("Rerender time:", performance.now() - renderTime);
      }
    });
    observer.observe(this.canvas);

    this.render();
  }
}



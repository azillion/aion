export class WebGPUCore {
  public readonly device: GPUDevice;
  public readonly context: GPUCanvasContext;
  public readonly presentationFormat: GPUTextureFormat;

  private constructor(device: GPUDevice, context: GPUCanvasContext, format: GPUTextureFormat) {
    this.device = device;
    this.context = context;
    this.presentationFormat = format;
  }

  static async create(canvas: HTMLCanvasElement): Promise<WebGPUCore> {
    if (!navigator.gpu) {
      throw new Error("WebGPU not supported on this browser.");
    }
    const adapter = await navigator.gpu.requestAdapter();
    if (!adapter) {
      throw new Error("No appropriate GPUAdapter found.");
    }

    const device = await adapter.requestDevice();
    const context = canvas.getContext("webgpu");
    if (!context) {
        throw new Error("Could not get WebGPU context from canvas.");
    }
    const presentationFormat = navigator.gpu.getPreferredCanvasFormat();

    context.configure({
      device,
      format: presentationFormat,
    });

    return new WebGPUCore(device, context, presentationFormat);
  }
}


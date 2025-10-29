import type { ISimulationHost } from '@shared/simulation';
import golShaderWGSL from '../../../game-sim/src/gol.wgsl?raw';

export class GpuSimHost implements ISimulationHost {
  private device!: GPUDevice;
  private queue!: GPUQueue;

  private stateA!: GPUTexture;
  private stateB!: GPUTexture;
  private pipeline!: GPUComputePipeline;
  private bindGroupAReadsBWrites!: GPUBindGroup;
  private bindGroupBReadsAWrites!: GPUBindGroup;
  private isStateACurrent = true;
  private gridSize = { width: 8192, height: 8192 };

  async initialize(device: GPUDevice, queue: GPUQueue): Promise<void> {
    this.device = device;
    this.queue = queue;

    const textureDescriptor: GPUTextureDescriptor = {
      size: this.gridSize,
      format: 'r32uint',
      usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.COPY_DST,
    };
    this.stateA = this.device.createTexture(textureDescriptor);
    this.stateB = this.device.createTexture(textureDescriptor);
    this.stateA.label = 'GOL State A';
    this.stateB.label = 'GOL State B';

    const initialData = new Uint32Array(this.gridSize.width * this.gridSize.height);
    for (let i = 0; i < initialData.length; i++) {
      initialData[i] = Math.random() < 0.25 ? 1 : 0;
    }
    this.queue.writeTexture(
      { texture: this.stateA },
      initialData,
      { bytesPerRow: this.gridSize.width * 4 },
      this.gridSize
    );

    const shaderModule = this.device.createShaderModule({ label: 'Game of Life Shader', code: golShaderWGSL });

    const bindGroupLayout = this.device.createBindGroupLayout({
      label: 'GOL Bind Group Layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'read-only', format: 'r32uint' } },
        { binding: 1, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'r32uint' } },
      ],
    });

    this.pipeline = await this.device.createComputePipelineAsync({
      label: 'GOL Compute Pipeline',
      layout: this.device.createPipelineLayout({ bindGroupLayouts: [bindGroupLayout] }),
      compute: { module: shaderModule, entryPoint: 'main' },
    });

    this.bindGroupAReadsBWrites = this.device.createBindGroup({
      label: 'GOL Bind Group (A -> B)',
      layout: bindGroupLayout,
      entries: [
        { binding: 0, resource: this.stateA.createView() },
        { binding: 1, resource: this.stateB.createView() },
      ],
    });
    this.bindGroupBReadsAWrites = this.device.createBindGroup({
      label: 'GOL Bind Group (B -> A)',
      layout: bindGroupLayout,
      entries: [
        { binding: 0, resource: this.stateB.createView() },
        { binding: 1, resource: this.stateA.createView() },
      ],
    });
  }

  tick(): void {
    const commandEncoder = this.device.createCommandEncoder();
    const pass = commandEncoder.beginComputePass();
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, this.isStateACurrent ? this.bindGroupAReadsBWrites : this.bindGroupBReadsAWrites);
    pass.dispatchWorkgroups(this.gridSize.width / 8, this.gridSize.height / 8);
    pass.end();
    this.queue.submit([commandEncoder.finish()]);
    this.isStateACurrent = !this.isStateACurrent;
  }

  getOutputTextureView(): GPUTextureView {
    return (this.isStateACurrent ? this.stateA : this.stateB).createView();
  }
}



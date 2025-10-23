import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import computeShaderWGSL from '../shaders/compute.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import { createShaderModule } from '../shaderUtils';

const WORKGROUP_SIZE = 8;

export class ComputePass implements IRenderPass {
  private pipeline!: GPUComputePipeline;
  private paramsUniformBuffer!: GPUBuffer;
  private core!: WebGPUCore; // Store core for re-initialization

  public async initialize(core: WebGPUCore, _scene: Scene): Promise<void> {
    this.core = core;
    const module = await createShaderModule(
      core.device,
      'Compute Shader Module',
      computeShaderWGSL,
      { 'camera.wgsl': cameraWGSL }
    );
    this.pipeline = this.core.device.createComputePipeline({
      label: 'Compute Pipeline',
      layout: 'auto',
      compute: {
        module,
        entryPoint: 'main'
      }
    });
    this.paramsUniformBuffer = core.device.createBuffer({
      size: 4,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }
  
  public async recreatePipeline(_bodyCount: number): Promise<void> {}

  public run(encoder: GPUCommandEncoder, context: RenderContext): void {
    const bodyCount = context.scene.lastKnownBodyCount;
    context.core.device.queue.writeBuffer(this.paramsUniformBuffer, 0, new Uint32Array([bodyCount]));

    const bindGroup = context.core.device.createBindGroup({
      label: "Compute Bind Group",
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: this.paramsUniformBuffer } },
        { binding: 1, resource: { buffer: context.scene.spheresBuffer } },
        { binding: 2, resource: context.mainSceneTexture.createView() },
        { binding: 3, resource: { buffer: context.scene.sharedCameraUniformBuffer } }
      ]
    });

    const pass = encoder.beginComputePass();
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    const workgroupCountX = Math.ceil(context.textureSize.width / WORKGROUP_SIZE);
    const workgroupCountY = Math.ceil(context.textureSize.height / WORKGROUP_SIZE);
    pass.dispatchWorkgroups(workgroupCountX, workgroupCountY);
    pass.end();
  }
}

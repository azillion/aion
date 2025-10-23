import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import computeShaderWGSL from '../shaders/compute.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';

const WORKGROUP_SIZE = 8;

export class ComputePass implements IRenderPass {
  private pipeline!: GPUComputePipeline;
  private core!: WebGPUCore; // Store core for re-initialization

  public initialize(core: WebGPUCore, scene: Scene): void {
    this.core = core;
    this.recreatePipeline(scene.lastKnownBodyCount);
  }
  
  public recreatePipeline(bodyCount: number): void {
    const module = this.core.device.createShaderModule({
      label: `Compute Shader Module (${bodyCount} bodies)`,
      code: `const NUM_SPHERES: u32 = ${bodyCount}u;\n` + cameraWGSL + computeShaderWGSL
    });

    this.pipeline = this.core.device.createComputePipeline({
      label: `Compute Pipeline (${bodyCount} bodies)`,
      layout: 'auto',
      compute: {
        module,
        entryPoint: 'main'
      }
    });
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext): void {
    const bindGroup = context.core.device.createBindGroup({
      label: "Compute Bind Group",
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: context.scene.spheresBuffer } },
        { binding: 1, resource: context.mainSceneTexture.createView() },
        { binding: 2, resource: { buffer: context.scene.sharedCameraUniformBuffer } }
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

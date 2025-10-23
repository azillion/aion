// imports kept for type parity; not used post-shared camera buffer
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import galaxyShaderWGSL from '../shaders/galaxy.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import { createShaderModule } from '../shaderUtils';

export class GalaxyPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private bindGroup!: GPUBindGroup;
  private depthTexture!: GPUTexture;

  public async initialize(core: WebGPUCore, scene: Scene): Promise<void> {
    const module = await createShaderModule(
      core.device,
      'Galaxy Shader Module',
      galaxyShaderWGSL,
      { 'camera.wgsl': cameraWGSL }
    );
    this.pipeline = core.device.createRenderPipeline({
      label: 'Galaxy Pipeline',
      layout: 'auto',
      vertex: { module, entryPoint: 'vertexMain' },
      fragment: { module, entryPoint: 'fragmentMain', targets: [{ format: core.presentationFormat }] },
      primitive: { topology: 'triangle-strip' },
      depthStencil: {
        format: 'depth24plus',
        depthWriteEnabled: true,
        depthCompare: 'less',
      },
    });

    this.bindGroup = core.device.createBindGroup({
        label: 'Galaxy Bind Group',
        layout: this.pipeline.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: { buffer: scene.sharedCameraUniformBuffer } },
          { binding: 1, resource: { buffer: scene.starBuffer } },
        ],
      });
  }
  
  public onResize(size: { width: number; height: number; }, core: WebGPUCore): void {
    if (this.depthTexture) this.depthTexture.destroy();
    this.depthTexture = core.device.createTexture({
        size,
        format: 'depth24plus',
        usage: GPUTextureUsage.RENDER_ATTACHMENT,
    });
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, starCount: number): void {
    // Shared camera buffer is already populated by Scene.update

    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: context.core.context.getCurrentTexture().createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: { r: 0.08, g: 0.0, b: 0.12, a: 1.0 }
      }],
      depthStencilAttachment: {
        view: this.depthTexture.createView(),
        depthClearValue: 1.0,
        depthLoadOp: 'clear',
        depthStoreOp: 'store',
      },
    });

    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, this.bindGroup);
    pass.draw(4, starCount);
    pass.end();
  }
}

import { mat4 } from 'gl-matrix';
import type { Vec3 } from '../../shared/types';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import galaxyShaderWGSL from '../shaders/galaxy.wgsl?raw';

export class GalaxyPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private bindGroup!: GPUBindGroup;
  private depthTexture!: GPUTexture;

  public initialize(core: WebGPUCore, scene: Scene): void {
    const module = core.device.createShaderModule({ code: galaxyShaderWGSL });
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
          { binding: 0, resource: { buffer: scene.galaxyCameraUniformBuffer } },
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
    // Update camera uniforms for this pass
    const proj = mat4.create();
    const aspect = context.textureSize.height > 0 ? (context.textureSize.width / context.textureSize.height) : 16/9;
    mat4.perspective(proj, Math.PI / 4, aspect, 0.1, 2000.0);
    const view = mat4.create();
    mat4.lookAt(view, context.camera.eye as any, context.camera.look_at as any, context.camera.up as any);
    const viewProj = mat4.multiply(mat4.create(), proj, view);
    const camRight: Vec3 = [view[0], view[4], view[8]];
    const camUp: Vec3 = [view[1], view[5], view[9]];
    const camData = new Float32Array(32);
    camData.set(viewProj, 0);
    camData.set(camRight, 16);
    camData.set(camUp, 20);
    context.core.device.queue.writeBuffer(context.scene.galaxyCameraUniformBuffer, 0, camData);

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

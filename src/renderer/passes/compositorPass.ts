import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import compositorWGSL from '../shaders/compositor.wgsl?raw';
import sceneUniformsWGSL from '../shaders/sceneUniforms.wgsl?raw';
import { createShaderModule } from '../shaderUtils';

export class CompositorPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private sampler!: GPUSampler;
  private bindGroupLayout!: GPUBindGroupLayout;

  public async initialize(core: WebGPUCore, _scene: Scene): Promise<void> {
    const module = await createShaderModule(core.device, 'Compositor Shader Module', compositorWGSL, {
      'sceneUniforms.wgsl': sceneUniformsWGSL,
    });
    this.bindGroupLayout = core.device.createBindGroupLayout({
      label: 'Compositor Bind Group Layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
        { binding: 1, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
        { binding: 2, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
        { binding: 3, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'unfilterable-float' } },
        { binding: 4, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'unfilterable-float' } },
        { binding: 5, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'unfilterable-float' } },
        { binding: 6, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
        { binding: 7, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
      ]
    });
    const layout = core.device.createPipelineLayout({ bindGroupLayouts: [this.bindGroupLayout] });
    this.pipeline = core.device.createRenderPipeline({
      label: 'Compositor Pipeline',
      layout,
      vertex: { module, entryPoint: 'vertexMain' },
      fragment: { module, entryPoint: 'fragmentMain', targets: [{ format: 'rgba16float' }] },
    });
    this.sampler = core.device.createSampler({ minFilter: 'linear', magFilter: 'linear' });
  }

  public run(
    encoder: GPUCommandEncoder,
    context: RenderContext,
    nearColor: GPUTexture,
    midColor: GPUTexture,
    farColor: GPUTexture,
    nearDepth: GPUTexture,
    midDepth: GPUTexture,
    farDepth: GPUTexture,
    destinationTexture: GPUTexture,
  ): void {
    const bindGroup = context.core.device.createBindGroup({
      label: 'Compositor Bind Group',
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: nearColor.createView() },
        { binding: 1, resource: midColor.createView() },
        { binding: 2, resource: farColor.createView() },
        { binding: 3, resource: nearDepth.createView() },
        { binding: 4, resource: midDepth.createView() },
        { binding: 5, resource: farDepth.createView() },
        { binding: 6, resource: this.sampler },
        { binding: 7, resource: { buffer: context.sceneUniformBuffer! } },
      ]
    });

    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: destinationTexture.createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: { r: 0, g: 0, b: 0, a: 1 },
      }]
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.draw(3);
    pass.end();
  }
}



import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import tierComputeWGSL from '../shaders/tierCompute.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import sceneUniformsWGSL from '../shaders/sceneUniforms.wgsl?raw';
import planetSdfWGSL from '../shaders/planetSdf.wgsl?raw';
import noiseWGSL from '../shaders/noise.wgsl?raw';
import atmosphereWGSL from '../shaders/atmosphere.wgsl?raw';
import { createShaderModule } from '../shaderUtils';

const WORKGROUP_SIZE = 8;

export class TierPass implements IRenderPass {
  private pipeline!: GPUComputePipeline;
  private paramsUniformBuffer!: GPUBuffer;
  private core!: WebGPUCore;
  private bindGroupLayout!: GPUBindGroupLayout;
  private tierBuffer!: GPUBuffer;
  private outputTexture!: GPUTexture;
  private depthTexture!: GPUTexture;

  constructor() {}

  public async initialize(core: WebGPUCore, _scene: Scene): Promise<void> {
    this.core = core;
    const module = await createShaderModule(
      core.device,
      'Tier Compute Shader Module',
      tierComputeWGSL,
      {
        'camera.wgsl': cameraWGSL,
        'sceneUniforms.wgsl': sceneUniformsWGSL,
        'planetSdf.wgsl': planetSdfWGSL,
        'noise.wgsl': noiseWGSL,
        'atmosphere.wgsl': atmosphereWGSL,
      }
    );
    this.bindGroupLayout = core.device.createBindGroupLayout({
      label: 'Tier Pass Bind Group Layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 2, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'rgba16float' } },
        { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 4, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'r32float' } },
        { binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 6, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 7, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
      ]
    });

    const pipelineLayout = core.device.createPipelineLayout({ bindGroupLayouts: [this.bindGroupLayout] });
    this.pipeline = this.core.device.createComputePipeline({
      label: 'Tier Compute Pipeline',
      layout: pipelineLayout,
      compute: { module, entryPoint: 'main' }
    });
    this.paramsUniformBuffer = core.device.createBuffer({ size: 4, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
  }

  public setOutputTexture(tex: GPUTexture) {
    this.outputTexture = tex;
  }
  public setDepthTexture(tex: GPUTexture) {
    this.depthTexture = tex;
  }
  public setTierBuffer(buf: GPUBuffer) {
    this.tierBuffer = buf;
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, bodyCount: number, sceneUniforms: GPUBuffer): void {
    if (!this.outputTexture || !this.depthTexture) return;
    context.core.device.queue.writeBuffer(this.paramsUniformBuffer, 0, new Uint32Array([bodyCount]));

    const bindGroup = context.core.device.createBindGroup({
      label: 'Tier Bind Group',
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: this.paramsUniformBuffer } },
        { binding: 1, resource: { buffer: this.tierBuffer } },
        { binding: 2, resource: this.outputTexture.createView() },
        { binding: 3, resource: { buffer: context.scene.sharedCameraUniformBuffer } },
        { binding: 4, resource: this.depthTexture.createView() },
        { binding: 5, resource: { buffer: sceneUniforms } },
        { binding: 6, resource: { buffer: context.scene.shadowCasterBuffer } },
        { binding: 7, resource: { buffer: context.scene.shadowCasterCountBuffer } },
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



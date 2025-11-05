import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import hydraulicsWGSL from './shaders/hydraulics.wgsl?raw';
import hydraulicsInitWGSL from './shaders/hydraulicsInit.wgsl?raw';
import cubeSphereWGSL from '../shaders/cubeSphere.wgsl?raw';
import { createShaderModule } from '../shaderUtils';

const WORKGROUP_SIZE = 8;
const SIM_GRID_SIZE = 1024; // Must match creation size in renderer

export class HydraulicsPass implements IRenderPass {
  private core!: WebGPUCore;
  private pipeline!: GPUComputePipeline;
  private bindGroupLayout!: GPUBindGroupLayout;
  private initPipeline!: GPUComputePipeline;

  public async initialize(core: WebGPUCore, _scene: Scene): Promise<void> {
    const module = await createShaderModule(
      core.device,
      'Hydraulics Compute Shader Module',
      hydraulicsWGSL,
      { 'cubeSphere.wgsl': cubeSphereWGSL }
    );
    const initModule = await createShaderModule(
      core.device,
      'Hydraulics Init Module',
      hydraulicsInitWGSL,
      { 'cubeSphere.wgsl': cubeSphereWGSL }
    );

    this.core = core;

    this.bindGroupLayout = core.device.createBindGroupLayout({
      label: 'Hydraulics BGL',
      entries: [
        { binding: 0, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'unfilterable-float', viewDimension: '2d-array' } },
        { binding: 1, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'rgba16float', viewDimension: '2d-array' } },
        { binding: 2, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'unfilterable-float', viewDimension: '2d-array' } },
      ]
    });

    const pipelineLayout = core.device.createPipelineLayout({
      label: 'Hydraulics Pipeline Layout',
      bindGroupLayouts: [this.bindGroupLayout]
    });

    this.pipeline = core.device.createComputePipeline({
      label: 'Hydraulics Compute Pipeline',
      layout: pipelineLayout,
      compute: { module, entryPoint: 'main' }
    });

    // Initialization pipeline (seed water based on terrain height and sea level)
    const initBGL = core.device.createBindGroupLayout({
      label: 'Hydraulics Init BGL',
      entries: [
        { binding: 0, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'rgba16float', viewDimension: '2d-array' } },
        { binding: 1, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'unfilterable-float', viewDimension: '2d-array' } },
        { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
      ]
    });
    const initPL = core.device.createPipelineLayout({ label: 'Hydraulics Init Pipeline Layout', bindGroupLayouts: [initBGL] });
    this.initPipeline = core.device.createComputePipeline({ label: 'Hydraulics Init Pipeline', layout: initPL, compute: { module: initModule, entryPoint: 'main' } });
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, terrainHeight: GPUTexture): void {
    const { core, waterRead, waterWrite } = context;
    if (!waterRead || !waterWrite) return;

    const bindGroup = core.device.createBindGroup({
      label: 'Hydraulics Bind Group',
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: waterRead.createView({ dimension: '2d-array' }) },
        { binding: 1, resource: waterWrite.createView({ dimension: '2d-array' }) },
        { binding: 2, resource: terrainHeight.createView({ dimension: '2d-array' }) },
      ]
    });

    const pass = encoder.beginComputePass();
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    const workgroupCountX = Math.ceil(SIM_GRID_SIZE / WORKGROUP_SIZE);
    const workgroupCountY = Math.ceil(SIM_GRID_SIZE / WORKGROUP_SIZE);
    pass.dispatchWorkgroups(workgroupCountX, workgroupCountY, 6);
    pass.end();
  }

  public initializeState(texture: GPUTexture, context: RenderContext, terrainHeight: GPUTexture): void {
    const initBGL = this.initPipeline.getBindGroupLayout(0);
    const encoder = this.core.device.createCommandEncoder();

    const planets = (context.scene as any).getPlanets ? (context.scene as any).getPlanets() : [];
    const planet = planets[0];
    if (!planet) return;
    const planetUBO = this.core.device.createBuffer({ size: 16, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    const u = new Float32Array(4);
    u[0] = planet.terrain!.radius;
    u[1] = planet.terrain!.seaLevel;
    this.core.device.queue.writeBuffer(planetUBO, 0, u);

    const bindGroup = this.core.device.createBindGroup({
      layout: initBGL,
      entries: [
        { binding: 0, resource: texture.createView({ dimension: '2d-array' }) },
        { binding: 1, resource: terrainHeight.createView({ dimension: '2d-array' }) },
        { binding: 2, resource: { buffer: planetUBO } },
      ]
    });

    const pass = encoder.beginComputePass();
    pass.setPipeline(this.initPipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(Math.ceil(texture.width / WORKGROUP_SIZE), Math.ceil(texture.height / WORKGROUP_SIZE), 6);
    pass.end();
    this.core.device.queue.submit([encoder.finish()]);

    setTimeout(() => planetUBO.destroy(), 100);
  }
}



import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import terrainBakeWGSL from './shaders/terrainBake.wgsl?raw';
import { createShaderModule } from '../shaderUtils';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import sceneUniformsWGSL from '../shaders/sceneUniforms.wgsl?raw';
import coarseGridMeshWGSL from '../shaders/coarseGridMesh.wgsl?raw';
import planetHeightWGSL from '../scene/shaders/planetHeight.wgsl?raw';
import noiseWGSL from '../shaders/noise.wgsl?raw';
import cubeSphereWGSL from '../shaders/cubeSphere.wgsl?raw';
import terrainCommonWGSL from '../scene/shaders/terrainCommon.wgsl?raw';

const WORKGROUP_SIZE = 8;

export class TerrainBakePass implements IRenderPass {
  private pipeline!: GPUComputePipeline;
  private bindGroupLayout!: GPUBindGroupLayout;
  private core!: WebGPUCore;

  public async initialize(core: WebGPUCore, _scene: Scene): Promise<void> {
    this.core = core;
    const module = await createShaderModule(core.device, 'Terrain Bake Module', terrainBakeWGSL, {
      'camera.wgsl': cameraWGSL,
      'sceneUniforms.wgsl': sceneUniformsWGSL,
      'coarseGrid.wgsl': coarseGridMeshWGSL,
      'planetHeight.wgsl': planetHeightWGSL,
      'noise.wgsl': noiseWGSL,
      'cubeSphere.wgsl': cubeSphereWGSL,
      'terrainCommon.wgsl': terrainCommonWGSL,
    });

    this.bindGroupLayout = core.device.createBindGroupLayout({
      label: 'Terrain Bake BGL',
      entries: [
        { binding: 0, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'r32float', viewDimension: '2d-array' } },
        { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 7, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 8, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 9, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
      ]
    });
    const pipelineLayout = core.device.createPipelineLayout({ bindGroupLayouts: [this.bindGroupLayout] });
    this.pipeline = core.device.createComputePipeline({
      label: 'Terrain Bake Pipeline',
      layout: pipelineLayout,
      compute: { module, entryPoint: 'main' },
    });
  }

  public run(context: RenderContext, outputTexture: GPUTexture): void {
    const { core, scene, gridVertexBuffer, gridElevationBuffer, gridIndexBuffer, sceneUniformBuffer } = context;
    const planets = (scene as any).getPlanets ? (scene as any).getPlanets() : [];
    const planet = planets[0];
    if (!planet || !gridVertexBuffer || !gridElevationBuffer || !sceneUniformBuffer) return;

    const paramsBuffer = core.device.createBuffer({ size: 16, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    const p = planet.terrain!;
    core.device.queue.writeBuffer(paramsBuffer, 0, new Float32Array([p.radius, p.seaLevel, p.maxHeight, p.noiseSeed]));

    const bindGroup = core.device.createBindGroup({
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: outputTexture.createView({ dimension: '2d-array' }) },
        { binding: 1, resource: { buffer: scene.sharedCameraUniformBuffer } },
        { binding: 2, resource: { buffer: sceneUniformBuffer } },
        { binding: 5, resource: { buffer: paramsBuffer } },
        { binding: 7, resource: { buffer: gridVertexBuffer } },
        { binding: 8, resource: { buffer: gridElevationBuffer } },
        { binding: 9, resource: { buffer: gridIndexBuffer! } },
      ]
    });

    const encoder = core.device.createCommandEncoder();
    const pass = encoder.beginComputePass();
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(Math.ceil(outputTexture.width / WORKGROUP_SIZE), Math.ceil(outputTexture.height / WORKGROUP_SIZE), 6);
    pass.end();
    core.device.queue.submit([encoder.finish()]);

    setTimeout(() => paramsBuffer.destroy(), 100);
  }
}



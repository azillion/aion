import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import sceneRenderWGSL from './shaders/sceneRender.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import cubeSphereWGSL from '../shaders/cubeSphere.wgsl?raw';
import sceneUniformsWGSL from '../shaders/sceneUniforms.wgsl?raw';
import coarseGridWGSL from '../shaders/coarseGrid.wgsl?raw';
import planetSdfWGSL from './shaders/planetSdf.wgsl?raw';
import noiseWGSL from '../shaders/noise.wgsl?raw';
import atmosphereWGSL from './shaders/atmosphere.wgsl?raw';
import raytracingWGSL from './shaders/raytracing.wgsl?raw';
import shadingWGSL from './shaders/shading.wgsl?raw';
import terrainCommonWGSL from './shaders/terrainCommon.wgsl?raw';
import f64WGSL from './shaders/f64.wgsl?raw';
import planetHeightWGSL from './shaders/planetHeight.wgsl?raw';
import { createShaderModule } from '../shaderUtils';

const WORKGROUP_SIZE = 8;

export class SceneRenderPass implements IRenderPass {
  private pipeline!: GPUComputePipeline;
  private bindGroupLayout!: GPUBindGroupLayout;
  private bodyCountBuffer!: GPUBuffer;

  public async initialize(core: WebGPUCore, _scene: Scene): Promise<void> {
    const module = await createShaderModule(
      core.device,
      'Scene Render Compute Shader Module',
      sceneRenderWGSL,
      {
        'camera.wgsl': cameraWGSL,
        'cubeSphere.wgsl': cubeSphereWGSL,
        'sceneUniforms.wgsl': sceneUniformsWGSL,
        'coarseGrid.wgsl': coarseGridWGSL,
        'terrainCommon.wgsl': terrainCommonWGSL,
        'planetSdf.wgsl': planetSdfWGSL,
        'planetHeight.wgsl': planetHeightWGSL,
        'noise.wgsl': noiseWGSL,
        'atmosphere.wgsl': atmosphereWGSL,
        'raytracing.wgsl': raytracingWGSL,
        'shading.wgsl': shadingWGSL,
        'f64.wgsl': f64WGSL,
      }
    );
    
    // Explicit bind group layout matching sceneRender.wgsl @group(0)
    this.bindGroupLayout = core.device.createBindGroupLayout({
      label: 'Scene Render BGL0',
      entries: [
        { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 2, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'rgba16float', viewDimension: '2d' } },
        { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 4, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 6, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 7, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 8, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 9, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 10, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'unfilterable-float', viewDimension: '2d-array' } },
      ]
    });

    const pipelineLayout = core.device.createPipelineLayout({
      label: 'Scene Render Pipeline Layout',
      bindGroupLayouts: [this.bindGroupLayout]
    });

    this.pipeline = core.device.createComputePipeline({
      label: 'Scene Render Compute Pipeline',
      layout: pipelineLayout,
      compute: { module, entryPoint: 'main' }
    });

    this.bodyCountBuffer = core.device.createBuffer({
      size: 4,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext): void {
    const { core, scene, mainSceneTexture } = context;

    core.device.queue.writeBuffer(this.bodyCountBuffer, 0, new Uint32Array([scene.sceneObjectCount]));

    const bindGroup = core.device.createBindGroup({
      label: 'Scene Render Bind Group',
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: this.bodyCountBuffer } },
        { binding: 1, resource: { buffer: scene.sceneObjectsBuffer } },
        { binding: 2, resource: mainSceneTexture.createView() },
        { binding: 3, resource: { buffer: scene.sharedCameraUniformBuffer } },
        { binding: 4, resource: { buffer: context.sceneUniformBuffer! } },
        { binding: 5, resource: { buffer: scene.shadowCasterBuffer } },
        { binding: 6, resource: { buffer: scene.shadowCasterCountBuffer } },
        { binding: 7, resource: { buffer: (context as any).gridVertexBuffer } },
        { binding: 8, resource: { buffer: (context as any).gridElevationBuffer } },
        { binding: 9, resource: { buffer: (context as any).gridIndexBuffer } },
        { binding: 10, resource: (context as any).waterWrite.createView({ dimension: '2d-array' }) },
      ]
    });

    const pass = encoder.beginComputePass();
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    const workgroupCountX = Math.ceil(context.textureSize.width / WORKGROUP_SIZE);
    const workgroupCountY = Math.ceil(context.textureSize.height / WORKGROUP_SIZE);
    pass.dispatchWorkgroups(workgroupCountX, workgroupCountY);
    pass.end();

    // No per-frame buffer destruction needed
  }
}



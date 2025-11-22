import type { RenderContext } from '../types';
import type { WebGPUCore } from '../core';
import gridOutlineWGSL from './shaders/gridOutline.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import { createShaderModule } from '../shaderUtils';

type PlanetInfo = {
  position: [number, number, number];
  radius: number;
};

export class GridOutlinePass {
  private pipeline!: GPURenderPipeline;
  private planetUniformBuffer!: GPUBuffer;
  private bindGroupLayout!: GPUBindGroupLayout;

  public async initialize(core: WebGPUCore): Promise<void> {
    const module = await createShaderModule(core.device, 'Grid Outline Shader', gridOutlineWGSL, {
      'camera.wgsl': cameraWGSL,
    });
    this.pipeline = core.device.createRenderPipeline({
      label: 'Grid Outline Pipeline',
      layout: 'auto',
      vertex: {
        module,
        entryPoint: 'vs_main',
        buffers: [
          {
            arrayStride: 12,
            attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x3' }],
          },
        ],
      },
      fragment: {
        module,
        entryPoint: 'fs_main',
        targets: [
          {
            format: 'rgba16float',
            blend: {
              color: { srcFactor: 'src-alpha', dstFactor: 'one-minus-src-alpha', operation: 'add' },
              alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add' },
            },
          },
        ],
      },
      primitive: {
        topology: 'line-list',
      },
    });
    this.bindGroupLayout = this.pipeline.getBindGroupLayout(0);
    this.planetUniformBuffer = core.device.createBuffer({
      label: 'Grid Outline Planet Uniform',
      size: 16,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }

  public run(
    encoder: GPUCommandEncoder,
    context: RenderContext,
    planet: PlanetInfo,
  ): void {
    if (
      !context.gridVertexBuffer ||
      !context.gridEdgeIndexBuffer ||
      !context.gridEdgeIndexCount ||
      context.gridEdgeIndexCount === 0
    ) {
      return;
    }
    const target = context.mainSceneTexture;
    if (!target) return;

    const device = context.core.device;
    const planetData = new Float32Array([
      planet.position[0],
      planet.position[1],
      planet.position[2],
      planet.radius,
    ]);
    device.queue.writeBuffer(this.planetUniformBuffer, 0, planetData);

    const bindGroup = device.createBindGroup({
      label: 'Grid Outline Bind Group',
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: context.scene.sharedCameraUniformBuffer } },
        { binding: 1, resource: { buffer: this.planetUniformBuffer } },
      ],
    });

    const pass = encoder.beginRenderPass({
      label: 'Grid Outline Pass',
      colorAttachments: [
        {
          view: target.createView(),
          loadOp: 'load',
          storeOp: 'store',
        },
      ],
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.setVertexBuffer(0, context.gridVertexBuffer);
    pass.setIndexBuffer(context.gridEdgeIndexBuffer, 'uint32');
    pass.drawIndexed(context.gridEdgeIndexCount);
    pass.end();
  }
}



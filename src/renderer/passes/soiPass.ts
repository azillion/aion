import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import soiShaderWGSL from '../shaders/soi.wgsl?raw';
import type { Body, Vec3 } from '../../shared/types';
import { G } from '../../shared/constants';

interface SoiInstance {
  position: Vec3;
  radius: number;
}

export class SOIPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private instanceBuffer!: GPUBuffer;
  private instanceData: SoiInstance[] = [];
  private instanceCapacity = 0;

  public initialize(core: WebGPUCore, _scene: Scene): void {
    const module = core.device.createShaderModule({ code: soiShaderWGSL });
    this.pipeline = core.device.createRenderPipeline({
      label: 'SOI Pipeline',
      layout: 'auto',
      vertex: { module, entryPoint: 'vertexMain' },
      fragment: {
        module,
        entryPoint: 'fragmentMain',
        targets: [{
          format: core.presentationFormat,
          blend: {
            color: { srcFactor: 'src-alpha', dstFactor: 'one-minus-src-alpha', operation: 'add' },
            alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add' },
          }
        }]
      },
      primitive: { topology: 'triangle-strip' },
    });
  }
  
  public update(bodies: Body[], scene: Scene, systemScale: number) {
    this.instanceData = [];
    const hierarchy = scene.hierarchy;

    for (const body of bodies) {
      const parentId = hierarchy.get(body.id);
      const parent = parentId ? bodies.find(b => b.id === parentId) : undefined;

      if (parent) {
        const r_vec = [body.position[0] - parent.position[0], body.position[1] - parent.position[1], body.position[2] - parent.position[2]] as Vec3;
        const v_vec = [body.velocity[0] - parent.velocity[0], body.velocity[1] - parent.velocity[1], body.velocity[2] - parent.velocity[2]] as Vec3;
        const r = Math.hypot(...r_vec);
        const v = Math.hypot(...v_vec);
        const mu = G * (parent.mass + body.mass);
        const invA = 2 / r - (v * v) / mu;
        if (invA <= 0) continue;
        const a = 1 / invA;

        const soiRadius = a * Math.pow(body.mass / parent.mass, 0.4);

        this.instanceData.push({
          position: [
            body.position[0] * systemScale,
            body.position[1] * systemScale,
            body.position[2] * systemScale,
          ],
          radius: soiRadius * systemScale,
        });
      }
    }
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext): void {
    if (this.instanceData.length === 0) return;
    const { core, scene } = context;

    if (this.instanceData.length > this.instanceCapacity) {
      if (this.instanceBuffer) this.instanceBuffer.destroy();
      this.instanceCapacity = Math.ceil(this.instanceData.length * 1.5);
      this.instanceBuffer = core.device.createBuffer({
        size: this.instanceCapacity * 16,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      });
    }

    const bufferData = new Float32Array(this.instanceData.length * 4);
    this.instanceData.forEach((inst, i) => {
      bufferData.set(inst.position, i * 4);
      bufferData[i * 4 + 3] = inst.radius;
    });
    core.device.queue.writeBuffer(this.instanceBuffer, 0, bufferData);

    const bindGroup = core.device.createBindGroup({
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: scene.mapCameraUniformBuffer } },
        { binding: 1, resource: { buffer: this.instanceBuffer } },
      ]
    });

    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: context.core.context.getCurrentTexture().createView(),
        loadOp: 'load',
        storeOp: 'store',
      }],
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.draw(4, this.instanceData.length);
    pass.end();
  }
}



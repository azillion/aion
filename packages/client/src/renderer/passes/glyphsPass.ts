import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import glyphsShaderWGSL from '../shaders/glyphs.wgsl?raw';
import type { Vec3 } from '@shared/types';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import { createShaderModule } from '../shaderUtils';

interface GlyphInstance {
  position: Vec3;
  scale: number;
}

export class GlyphsPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private instanceBuffer!: GPUBuffer;
  private instanceData: GlyphInstance[] = [];
  private instanceCapacity = 0;

  public async initialize(core: WebGPUCore, _scene: Scene): Promise<void> {
    const module = await createShaderModule(
      core.device,
      'Glyphs Shader Module',
      glyphsShaderWGSL,
      { 'camera.wgsl': cameraWGSL }
    );
    this.pipeline = core.device.createRenderPipeline({
      label: 'Glyphs Pipeline',
      layout: 'auto',
      vertex: { module, entryPoint: 'vertexMain' },
      fragment: { module, entryPoint: 'fragmentMain', targets: [{ format: core.presentationFormat }] },
    });
  }
  
  public update(periapsisPoints: (Vec3 | null)[], apoapsisPoints: (Vec3 | null)[], ascendingNodePoints: (Vec3 | null)[], descendingNodePoints: (Vec3 | null)[]) {
    this.instanceData = [];
    for (const p of periapsisPoints) {
        if (!p) continue;
        this.instanceData.push({ position: p, scale: 6.0 }); // upright
    }
    for (const p of apoapsisPoints) {
        if (!p) continue;
        this.instanceData.push({ position: p, scale: -6.0 }); // inverted
    }
    for (const p of ascendingNodePoints) {
        if (!p) continue;
        this.instanceData.push({ position: p, scale: 6.0 });
    }
    for (const p of descendingNodePoints) {
        if (!p) continue;
        this.instanceData.push({ position: p, scale: -6.0 });
    }
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext): void {
    if (this.instanceData.length === 0) return;
    const { core, scene } = context;

    // Resize buffer if needed
    if (this.instanceData.length > this.instanceCapacity) {
        if (this.instanceBuffer) this.instanceBuffer.destroy();
        this.instanceCapacity = Math.ceil(this.instanceData.length * 1.5);
        this.instanceBuffer = core.device.createBuffer({
            size: this.instanceCapacity * 16, // 4 floats * 4 bytes (vec3 + scale)
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
        });
    }

    // Write instance data
    const bufferData = new Float32Array(this.instanceData.length * 4);
    this.instanceData.forEach((inst, i) => {
        const base = i * 4;
        bufferData.set(inst.position, base);
        bufferData[base + 3] = inst.scale;
    });
    core.device.queue.writeBuffer(this.instanceBuffer, 0, bufferData);

    const bindGroup = core.device.createBindGroup({
        layout: this.pipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: scene.sharedCameraUniformBuffer } },
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
    pass.draw(3, this.instanceData.length);
    pass.end();
  }
}



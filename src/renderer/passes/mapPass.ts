import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import mapShaderWGSL from '../shaders/map.wgsl?raw';

export class MapPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;

  public initialize(core: WebGPUCore, _scene: Scene): void {
    const module = core.device.createShaderModule({ code: mapShaderWGSL });
    this.pipeline = core.device.createRenderPipeline({
      label: 'Map Pipeline',
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

  public run(encoder: GPUCommandEncoder, context: RenderContext): void {
    const { core, scene, camera } = context;

    // 1. Get the pre-computed view-projection matrix from the camera
    const viewProj = ((): Float32Array => {
      const m = (window as any)._tmp_m || new Float32Array(16);
      (window as any)._tmp_m = m;
      // m = projection * view
      // manual multiply to avoid importing gl-matrix here
      const a = camera.projectionMatrix as unknown as Float32Array;
      const b = camera.viewMatrix as unknown as Float32Array;
      for (let i = 0; i < 4; i++) {
        const ai0 = a[i*4+0], ai1 = a[i*4+1], ai2 = a[i*4+2], ai3 = a[i*4+3];
        m[i*4+0] = ai0*b[0] + ai1*b[4] + ai2*b[8]  + ai3*b[12];
        m[i*4+1] = ai0*b[1] + ai1*b[5] + ai2*b[9]  + ai3*b[13];
        m[i*4+2] = ai0*b[2] + ai1*b[6] + ai2*b[10] + ai3*b[14];
        m[i*4+3] = ai0*b[3] + ai1*b[7] + ai2*b[11] + ai3*b[15];
      }
      return m;
    })();
    core.device.queue.writeBuffer(scene.mapCameraUniformBuffer, 0, viewProj as unknown as ArrayBuffer);

    // 2. Create Bind Group
    const bindGroup = core.device.createBindGroup({
      label: 'Map Pass Bind Group',
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: scene.mapCameraUniformBuffer } },
        { binding: 1, resource: { buffer: scene.spheresBuffer } },
      ],
    });

    // 3. Render
    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: core.context.getCurrentTexture().createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: { r: 0.01, g: 0.02, b: 0.03, a: 1.0 }
      }],
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.draw(4, scene.lastKnownBodyCount);
    pass.end();
  }
}



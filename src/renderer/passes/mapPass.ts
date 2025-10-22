import { mat4 } from 'gl-matrix';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import type { Theme } from '../../shared/types';
import mapShaderWGSL from '../shaders/map.wgsl?raw';
import backgroundShaderWGSL from '../shaders/background.wgsl?raw';

export class MapPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private backgroundPipeline!: GPURenderPipeline;

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

    // Background pipeline to draw the precomputed orbits texture
    const bgModule = core.device.createShaderModule({ code: backgroundShaderWGSL });
    this.backgroundPipeline = core.device.createRenderPipeline({
        label: 'Map Background Pipeline',
        layout: 'auto',
        vertex: { module: bgModule, entryPoint: 'vertexMain' },
        fragment: {
            module: bgModule,
            entryPoint: 'fragmentMain',
            targets: [{ format: core.presentationFormat }],
        },
    });
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, _theme: Theme, showOrbits: boolean): void {
    const { core, scene, camera, orbitsTexture } = context;

    // 1. Build camera uniform (viewProjection, right, up)
    const viewProj = mat4.multiply(mat4.create(), camera.projectionMatrix, camera.viewMatrix);
    const eye = camera.eye as unknown as number[];
    const look = camera.look_at as unknown as number[];
    const upWorld = camera.up as unknown as number[];
    const fwd = [
      look[0] - eye[0],
      look[1] - eye[1],
      look[2] - eye[2],
    ];
    const fwdLen = Math.hypot(fwd[0], fwd[1], fwd[2]) || 1.0;
    const f = [fwd[0]/fwdLen, fwd[1]/fwdLen, fwd[2]/fwdLen];
    // right = normalize(cross(fwd, up))
    const r = [
      f[1]*upWorld[2] - f[2]*upWorld[1],
      f[2]*upWorld[0] - f[0]*upWorld[2],
      f[0]*upWorld[1] - f[1]*upWorld[0],
    ];
    const rLen = Math.hypot(r[0], r[1], r[2]) || 1.0;
    const right = [r[0]/rLen, r[1]/rLen, r[2]/rLen];
    // recompute up to ensure orthonormality: up = normalize(cross(right, fwd))
    const u = [
      right[1]*f[2] - right[2]*f[1],
      right[2]*f[0] - right[0]*f[2],
      right[0]*f[1] - right[1]*f[0],
    ];
    const uLen = Math.hypot(u[0], u[1], u[2]) || 1.0;
    const up = [u[0]/uLen, u[1]/uLen, u[2]/uLen];

    const camData = new Float32Array(28);
    camData.set(viewProj, 0);
    camData.set([right[0], right[1], right[2], 0.0], 16);
    camData.set([up[0], up[1], up[2], 0.0], 20);
    // last vec4 padding remains zero
    core.device.queue.writeBuffer(scene.mapCameraUniformBuffer, 0, camData);

    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: core.context.getCurrentTexture().createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: { r: 0.0, g: 0.0, b: 0.0, a: 1.0 }
      }],
    });
    
    // 2. Draw the Orbits Texture as the background
    if (showOrbits) {
        const sampler = core.device.createSampler();
        const backgroundBindGroup = core.device.createBindGroup({
            label: 'Map Background BindGroup',
            layout: this.backgroundPipeline.getBindGroupLayout(0),
            entries: [
                { binding: 0, resource: sampler },
                { binding: 1, resource: orbitsTexture.createView() }
            ]
        });
        pass.setPipeline(this.backgroundPipeline);
        pass.setBindGroup(0, backgroundBindGroup);
        pass.draw(6);
    }

    // 3. Draw the planet icons on top
    const planetBindGroup = core.device.createBindGroup({
      label: 'Map Planet Bind Group',
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: scene.mapCameraUniformBuffer } },
        { binding: 1, resource: { buffer: scene.spheresBuffer } },
      ],
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, planetBindGroup);
    pass.draw(4, scene.lastKnownBodyCount);
    pass.end();
  }
}



import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import postfxShaderWGSL from '../shaders/postfx.wgsl?raw';
import type { Body, Theme } from '../../shared/types';

export class PostFXPass implements IRenderPass {
  private postfxPipeline!: GPURenderPipeline;
  private presentPipeline!: GPURenderPipeline;
  private sampler!: GPUSampler;
  private targetInfoUniform!: GPUBuffer;

  public initialize(core: WebGPUCore, _scene: Scene): void {
    const module = core.device.createShaderModule({ code: postfxShaderWGSL });

    this.postfxPipeline = core.device.createRenderPipeline({
      label: 'PostFX Pipeline',
      layout: "auto",
      vertex: { module, entryPoint: "vertexMain" },
      fragment: { module, entryPoint: "fragmentMain", targets: [{ format: core.presentationFormat }] },
    });
    
    this.presentPipeline = core.device.createRenderPipeline({
        label: 'Present Pipeline',
        layout: "auto",
        vertex: { module, entryPoint: "vertexMain" },
        fragment: { module, entryPoint: "presentFragment", targets: [{ format: core.presentationFormat }] },
    });

    this.sampler = core.device.createSampler({ magFilter: "linear", minFilter: "linear" });

    this.targetInfoUniform = core.device.createBuffer({
      label: 'Target Info Uniform',
      size: 272, // Matches WGSL OrbitMask struct
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, theme: Theme, bodies: Body[]): void {
    const { core, destinationTexture, sourceTexture } = context;

    this.updateTargetInfo(context, bodies);

    const postFxBindGroup = core.device.createBindGroup({
      label: 'PostFX Bind Group',
      layout: this.postfxPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: this.sampler },
        { binding: 1, resource: context.mainSceneTexture.createView() },
        { binding: 2, resource: { buffer: context.themeUniformBuffer } },
        { binding: 3, resource: sourceTexture.createView() },
        { binding: 4, resource: context.orbitsTexture.createView() },
        { binding: 5, resource: { buffer: this.targetInfoUniform } },
      ]
    });

    // Pass 1: Render PostFX into offscreen destination texture
    const postPass = encoder.beginRenderPass({
      colorAttachments: [{
        view: destinationTexture.createView(),
        loadOp: "clear",
        storeOp: "store",
        clearValue: { r: 0, g: 0, b: 0, a: 1 }
      }]
    });
    postPass.setPipeline(this.postfxPipeline);
    postPass.setBindGroup(0, postFxBindGroup);
    postPass.draw(6);
    postPass.end();

    // Pass 2: Present the destination texture to the canvas
    const presentBindGroup = core.device.createBindGroup({
      label: 'Present Bind Group',
      layout: this.presentPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: this.sampler },
        { binding: 1, resource: destinationTexture.createView() },
      ]
    });
    const presentPass = encoder.beginRenderPass({
      colorAttachments: [{
        view: core.context.getCurrentTexture().createView(),
        loadOp: "clear",
        storeOp: "store",
        clearValue: { r: 0, g: 0, b: 0, a: 1 }
      }]
    });
    presentPass.setPipeline(this.presentPipeline);
    presentPass.setBindGroup(0, presentBindGroup);
    presentPass.draw(6);
    presentPass.end();
  }
  
  private updateTargetInfo(context: RenderContext, bodies: Body[]) {
    const { core, camera, systemScale, textureSize } = context;
    const buf = new ArrayBuffer(272);
    const u32 = new Uint32Array(buf);
    const f32 = new Float32Array(buf);
    
    u32[0] = bodies.length;
    f32[1] = 20.0; // targetRadiusPx
    f32[2] = textureSize.width;
    f32[3] = textureSize.height;

    for (let i=0; i < bodies.length; i++) {
        const b = bodies[i];
        const worldPos = [ b.position[0] * systemScale, b.position[1] * systemScale, b.position[2] * systemScale ];
        const camRelative = [worldPos[0] - camera.eye[0], worldPos[1] - camera.eye[1], worldPos[2] - camera.eye[2]];
        const lookDir = [camera.look_at[0] - camera.eye[0], camera.look_at[1] - camera.eye[1], camera.look_at[2] - camera.eye[2]];
        const lookLen = Math.hypot(...lookDir) || 1;
        const forward = [lookDir[0]/lookLen, lookDir[1]/lookLen, lookDir[2]/lookLen];
        const up = camera.up;
        const r_x = up[1]*forward[2] - up[2]*forward[1], r_y = up[2]*forward[0] - up[0]*forward[2], r_z = up[0]*forward[1] - up[1]*forward[0];
        const rightLen = Math.hypot(r_x, r_y, r_z) || 1;
        const right = [r_x/rightLen, r_y/rightLen, r_z/rightLen];
        const camUp = [forward[1]*right[2] - forward[2]*right[1], forward[2]*right[0] - forward[0]*right[2], forward[0]*right[1] - forward[1]*right[0]];
        const viewX = camRelative[0]*right[0] + camRelative[1]*right[1] + camRelative[2]*right[2];
        const viewY = camRelative[0]*camUp[0] + camRelative[1]*camUp[1] + camRelative[2]*camUp[2];
        const viewZ = camRelative[0]*forward[0] + camRelative[1]*forward[1] + camRelative[2]*forward[2];

        if (viewZ <= 0) continue;

        const vfovRad = 25 * Math.PI / 180;
        const projScale = 1 / Math.tan(vfovRad / 2.0);
        const ndcX = (viewX / viewZ) * projScale;
        const ndcY = (viewY / viewZ) * projScale;
        const px = (ndcX * 0.5 + 0.5) * textureSize.width;
        const py = (1.0 - (ndcY * 0.5 + 0.5)) * textureSize.height;

        const base = 4 + i * 4;
        f32[base + 0] = px;
        f32[base + 1] = py;
    }
    
    core.device.queue.writeBuffer(this.targetInfoUniform, 0, buf);
  }

  public static clearTexture(device: GPUDevice, texture: GPUTexture, color = { r: 0, g: 0, b: 0, a: 1 }) {
    const encoder = device.createCommandEncoder();
    encoder.beginRenderPass({
      colorAttachments: [{ view: texture.createView(), loadOp: 'clear', storeOp: 'store', clearValue: color }]
    }).end();
    device.queue.submit([encoder.finish()]);
  }
}

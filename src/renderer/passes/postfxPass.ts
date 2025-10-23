import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import postfxShaderWGSL from '../shaders/postfx.wgsl?raw';
import type { Body, Theme, Vec3 } from '../../shared/types';
import { projectWorldToScreen } from '../projection';

export class PostFXPass implements IRenderPass {
  private postfxPipeline!: GPURenderPipeline;
  private presentPipeline!: GPURenderPipeline;
  private sampler!: GPUSampler;
  private targetInfoUniform!: GPUBuffer;
  private bindGroupLayoutPost!: GPUBindGroupLayout;
  private bindGroupLayoutPresent!: GPUBindGroupLayout;

  public initialize(core: WebGPUCore, _scene: Scene): void {
    const module = core.device.createShaderModule({ code: postfxShaderWGSL });

    // Explicit layouts: one for full post-fx, one for simple present
    this.bindGroupLayoutPost = core.device.createBindGroupLayout({
      label: 'PostFX Bind Group Layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
        { binding: 1, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
        { binding: 2, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
        { binding: 3, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
        { binding: 4, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
        { binding: 5, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
      ]
    });
    const postPipelineLayout = core.device.createPipelineLayout({ bindGroupLayouts: [this.bindGroupLayoutPost] });

    this.bindGroupLayoutPresent = core.device.createBindGroupLayout({
      label: 'Present Bind Group Layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
        { binding: 1, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
      ]
    });
    const presentPipelineLayout = core.device.createPipelineLayout({ bindGroupLayouts: [this.bindGroupLayoutPresent] });

    this.postfxPipeline = core.device.createRenderPipeline({
      label: 'PostFX Pipeline',
      layout: postPipelineLayout,
      vertex: { module, entryPoint: "vertexMain" },
      fragment: { module, entryPoint: "fragmentMain", targets: [{ format: core.presentationFormat }] },
    });
    
    this.presentPipeline = core.device.createRenderPipeline({
        label: 'Present Pipeline',
        layout: presentPipelineLayout,
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

  public getPresentPipeline(): GPURenderPipeline {
    return this.presentPipeline;
  }

  public getSampler(): GPUSampler {
    return this.sampler;
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, theme: Theme, bodies: Body[]): void {
    const { core, destinationTexture, sourceTexture } = context;

    this.updateTargetInfo(context, bodies);

    const postFxBindGroup = core.device.createBindGroup({
      label: 'PostFX Bind Group',
      layout: this.bindGroupLayoutPost,
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
      layout: this.bindGroupLayoutPresent,
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

    for (let i = 0; i < bodies.length; i++) {
        const b = bodies[i];
        const worldPos: Vec3 = [
            b.position[0] * systemScale,
            b.position[1] * systemScale,
            b.position[2] * systemScale
        ];

        const screenPos = projectWorldToScreen(worldPos, camera, textureSize);

        if (screenPos) {
            const base = 4 + i * 4;
            f32[base + 0] = screenPos.x;
            f32[base + 1] = screenPos.y;
        }
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

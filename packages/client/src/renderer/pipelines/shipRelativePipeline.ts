import type { FrameData, Theme } from '@shared/types';
import type { RenderContext } from '../types';
import type { IRenderPipeline } from './base';
import { FAR_TIER_SCALE, MID_TIER_SCALE } from '@shared/constants';
import { PostFXPass } from '../passes/postfxPass';
import { TierPass } from '../passes/tierPass';
import { CompositorPass } from '../passes/compositorPass';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';

export class ShipRelativePipeline implements IRenderPipeline {
  private nearTierPass: TierPass;
  private midTierPass: TierPass;
  private farTierPass: TierPass;
  private compositorPass: CompositorPass;
  private postfxPass: PostFXPass;

  constructor() {
    this.nearTierPass = new TierPass();
    this.midTierPass = new TierPass();
    this.farTierPass = new TierPass();
    this.compositorPass = new CompositorPass();
    this.postfxPass = new PostFXPass();
  }

  private nearColorTexture!: GPUTexture;
  private nearDepthTexture!: GPUTexture;
  private midColorTexture!: GPUTexture;
  private midDepthTexture!: GPUTexture;
  private farColorTexture!: GPUTexture;
  private farDepthTexture!: GPUTexture;
  private nearSceneUniformBuffer!: GPUBuffer;
  private midSceneUniformBuffer!: GPUBuffer;
  private farSceneUniformBuffer!: GPUBuffer;

  public async initialize(core: WebGPUCore, scene: Scene): Promise<void> {
    const sceneUniformBufferSize = 48; // three vec4<f32>
    this.nearSceneUniformBuffer = core.device.createBuffer({
      size: sceneUniformBufferSize,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.midSceneUniformBuffer = core.device.createBuffer({
      size: sceneUniformBufferSize,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.farSceneUniformBuffer = core.device.createBuffer({
      size: sceneUniformBufferSize,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    await this.nearTierPass.initialize(core, scene);
    await this.midTierPass.initialize(core, scene);
    await this.farTierPass.initialize(core, scene);
    await this.compositorPass.initialize(core, scene);
    await this.postfxPass.initialize(core, scene);
  }

  public onResize(size: { width: number; height: number }, core: WebGPUCore): void {
    [
      this.nearColorTexture, this.nearDepthTexture,
      this.midColorTexture, this.midDepthTexture,
      this.farColorTexture, this.farDepthTexture
    ].forEach(tex => tex?.destroy());

    const colorUsage = GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.STORAGE_BINDING;
    this.nearColorTexture = core.device.createTexture({ size, format: 'rgba16float', usage: colorUsage });
    this.nearDepthTexture = core.device.createTexture({ size, format: 'r32float', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING });
    this.midColorTexture = core.device.createTexture({ size, format: 'rgba16float', usage: colorUsage });
    this.midDepthTexture = core.device.createTexture({ size, format: 'r32float', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING });
    this.farColorTexture = core.device.createTexture({ size, format: 'rgba16float', usage: colorUsage });
    this.farDepthTexture = core.device.createTexture({ size, format: 'r32float', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING });
  }

  render(encoder: GPUCommandEncoder, context: RenderContext, frameData: FrameData, theme: Theme): void {
    const clearBlack = { r: 0, g: 0, b: 0, a: 1 };
    encoder.beginRenderPass({
      colorAttachments: [{ view: this.nearColorTexture.createView(), loadOp: 'clear', storeOp: 'store', clearValue: clearBlack }]
    }).end();
    encoder.beginRenderPass({
      colorAttachments: [{ view: this.midColorTexture.createView(), loadOp: 'clear', storeOp: 'store', clearValue: clearBlack }]
    }).end();
    encoder.beginRenderPass({
      colorAttachments: [{ view: this.farColorTexture.createView(), loadOp: 'clear', storeOp: 'store', clearValue: clearBlack }]
    }).end();

    this.nearTierPass.setOutputTexture(this.nearColorTexture);
    this.nearTierPass.setDepthTexture(this.nearDepthTexture);
    this.midTierPass.setOutputTexture(this.midColorTexture);
    this.midTierPass.setDepthTexture(this.midDepthTexture);
    this.farTierPass.setOutputTexture(this.farColorTexture);
    this.farTierPass.setDepthTexture(this.farDepthTexture);
    this.nearTierPass.setTierBuffer(context.scene.nearTierBuffer);
    this.midTierPass.setTierBuffer(context.scene.midTierBuffer);
    this.farTierPass.setTierBuffer(context.scene.farTierBuffer);

    if (frameData.dominantLight && frameData.worldCameraEye) {
      const light = frameData.dominantLight;
      const emissive: [number, number, number] = [1.0, 1.0, 1.0];
      const LIGHT_INTENSITY = 15.0;

      let referenceBody = frameData.rawState.bodies.find(b => b.id === 'sol');
      let closestDistSq = Number.POSITIVE_INFINITY;
      for (const b of frameData.rawState.bodies) {
        if (typeof b.mass === 'number' && b.mass >= 1e22) {
          const dx = b.position[0] - frameData.worldCameraEye[0];
          const dy = b.position[1] - frameData.worldCameraEye[1];
          const dz = b.position[2] - frameData.worldCameraEye[2];
          const d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < closestDistSq) { closestDistSq = d2; referenceBody = b; }
        }
      }
      const refPos = referenceBody ? referenceBody.position : ([0,0,0] as [number, number, number]);
      const lightDir: [number, number, number] = [
        refPos[0] - light.position[0],
        refPos[1] - light.position[1],
        refPos[2] - light.position[2],
      ];
      const dirLen = Math.hypot(lightDir[0], lightDir[1], lightDir[2]);
      if (dirLen > 0) { lightDir[0] /= dirLen; lightDir[1] /= dirLen; lightDir[2] /= dirLen; }

      this.updateTierLightingUniforms(context.core, frameData, this.farSceneUniformBuffer, lightDir, emissive, LIGHT_INTENSITY, frameData.debugTierView ?? -1, FAR_TIER_SCALE);
      this.farTierPass.run(encoder, context, context.scene.farCount, this.farSceneUniformBuffer);
      this.updateTierLightingUniforms(context.core, frameData, this.midSceneUniformBuffer, lightDir, emissive, LIGHT_INTENSITY, frameData.debugTierView ?? -1, MID_TIER_SCALE);
      this.midTierPass.run(encoder, context, context.scene.midCount, this.midSceneUniformBuffer);
      this.updateTierLightingUniforms(context.core, frameData, this.nearSceneUniformBuffer, lightDir, emissive, LIGHT_INTENSITY, frameData.debugTierView ?? -1, 1.0);
      this.nearTierPass.run(encoder, context, context.scene.nearCount, this.nearSceneUniformBuffer);
    } else {
      this.updateTierLightingUniforms(context.core, frameData, this.farSceneUniformBuffer, [0,0,-1], [0,0,0], 0, frameData.debugTierView ?? -1, FAR_TIER_SCALE);
      this.farTierPass.run(encoder, context, context.scene.farCount, this.farSceneUniformBuffer);
      this.updateTierLightingUniforms(context.core, frameData, this.midSceneUniformBuffer, [0,0,-1], [0,0,0], 0, frameData.debugTierView ?? -1, MID_TIER_SCALE);
      this.midTierPass.run(encoder, context, context.scene.midCount, this.midSceneUniformBuffer);
      this.updateTierLightingUniforms(context.core, frameData, this.nearSceneUniformBuffer, [0,0,-1], [0,0,0], 0, frameData.debugTierView ?? -1, 1.0);
      this.nearTierPass.run(encoder, context, context.scene.nearCount, this.nearSceneUniformBuffer);
    }

    this.compositorPass.run(
      encoder,
      context,
      this.nearColorTexture,
      this.midColorTexture,
      this.farColorTexture,
      this.nearDepthTexture,
      this.midDepthTexture,
      this.farDepthTexture,
      context.mainSceneTexture,
    );

    PostFXPass.clearTexture(context.core.device, context.orbitsTexture);
    this.postfxPass.run(encoder, context, theme, frameData.rawState.bodies ?? []);
  }

  private updateTierLightingUniforms(
    core: WebGPUCore,
    frameData: FrameData,
    targetBuffer: GPUBuffer,
    lightDirection: [number, number, number],
    emissive: [number, number, number],
    intensity: number,
    debugTierView: number,
    tierScale: number
  ) {
    const len = Math.hypot(emissive[0], emissive[1], emissive[2]);
    const finalLightColor: [number, number, number] = len > 0
      ? [
          (emissive[0] / len) * intensity,
          (emissive[1] / len) * intensity,
          (emissive[2] / len) * intensity,
        ]
      : [0, 0, 0];
    const bufferData = new Float32Array(12);
    bufferData.set(lightDirection, 0);
    bufferData[3] = 0.0;
    bufferData.set(finalLightColor, 4);
    bufferData[7] = debugTierView;
    bufferData[8] = tierScale;
    bufferData[9] = frameData.showAtmosphere ? 1.0 : 0.0;
    bufferData[10] = 0.0;
    bufferData[11] = 0.0;
    core.device.queue.writeBuffer(targetBuffer, 0, bufferData);
  }

  public prepare(_frameData: FrameData, _renderer: any): void {
    // No-op for this pipeline
  }
}



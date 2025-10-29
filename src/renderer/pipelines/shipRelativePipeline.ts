import type { FrameData, Theme } from '../../shared/types';
import type { RenderContext } from '../types';
import type { IRenderPipeline } from './base';
import { FAR_TIER_SCALE, MID_TIER_SCALE } from '../../shared/constants';
import { PostFXPass } from '../passes/postfxPass';

export class ShipRelativePipeline implements IRenderPipeline {
  constructor(private renderer: any) {}

  render(encoder: GPUCommandEncoder, context: RenderContext, frameData: FrameData, theme: Theme): void {
    const clearBlack = { r: 0, g: 0, b: 0, a: 1 };
    encoder.beginRenderPass({
      colorAttachments: [{ view: this.renderer.nearColorTexture.createView(), loadOp: 'clear', storeOp: 'store', clearValue: clearBlack }]
    }).end();
    encoder.beginRenderPass({
      colorAttachments: [{ view: this.renderer.midColorTexture.createView(), loadOp: 'clear', storeOp: 'store', clearValue: clearBlack }]
    }).end();
    encoder.beginRenderPass({
      colorAttachments: [{ view: this.renderer.farColorTexture.createView(), loadOp: 'clear', storeOp: 'store', clearValue: clearBlack }]
    }).end();

    this.renderer.nearTierPass.setOutputTexture(this.renderer.nearColorTexture);
    this.renderer.nearTierPass.setDepthTexture(this.renderer.nearDepthTexture);
    this.renderer.midTierPass.setOutputTexture(this.renderer.midColorTexture);
    this.renderer.midTierPass.setDepthTexture(this.renderer.midDepthTexture);
    this.renderer.farTierPass.setOutputTexture(this.renderer.farColorTexture);
    this.renderer.farTierPass.setDepthTexture(this.renderer.farDepthTexture);
    this.renderer.nearTierPass.setTierBuffer(this.renderer.scene.nearTierBuffer);
    this.renderer.midTierPass.setTierBuffer(this.renderer.scene.midTierBuffer);
    this.renderer.farTierPass.setTierBuffer(this.renderer.scene.farTierBuffer);

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

      this.renderer.updateTierLightingUniforms(this.renderer.farSceneUniformBuffer, lightDir, emissive, LIGHT_INTENSITY, frameData.debugTierView ?? -1, FAR_TIER_SCALE);
      this.renderer.farTierPass.run(encoder, context, this.renderer.scene.farCount, this.renderer.farSceneUniformBuffer);
      this.renderer.updateTierLightingUniforms(this.renderer.midSceneUniformBuffer, lightDir, emissive, LIGHT_INTENSITY, frameData.debugTierView ?? -1, MID_TIER_SCALE);
      this.renderer.midTierPass.run(encoder, context, this.renderer.scene.midCount, this.renderer.midSceneUniformBuffer);
      this.renderer.updateTierLightingUniforms(this.renderer.nearSceneUniformBuffer, lightDir, emissive, LIGHT_INTENSITY, frameData.debugTierView ?? -1, 1.0);
      this.renderer.nearTierPass.run(encoder, context, this.renderer.scene.nearCount, this.renderer.nearSceneUniformBuffer);
    } else {
      this.renderer.updateTierLightingUniforms(this.renderer.farSceneUniformBuffer, [0,0,-1], [0,0,0], 0, frameData.debugTierView ?? -1, FAR_TIER_SCALE);
      this.renderer.farTierPass.run(encoder, context, this.renderer.scene.farCount, this.renderer.farSceneUniformBuffer);
      this.renderer.updateTierLightingUniforms(this.renderer.midSceneUniformBuffer, [0,0,-1], [0,0,0], 0, frameData.debugTierView ?? -1, MID_TIER_SCALE);
      this.renderer.midTierPass.run(encoder, context, this.renderer.scene.midCount, this.renderer.midSceneUniformBuffer);
      this.renderer.updateTierLightingUniforms(this.renderer.nearSceneUniformBuffer, [0,0,-1], [0,0,0], 0, frameData.debugTierView ?? -1, 1.0);
      this.renderer.nearTierPass.run(encoder, context, this.renderer.scene.nearCount, this.renderer.nearSceneUniformBuffer);
    }

    this.renderer.compositorPass.run(
      encoder,
      context,
      this.renderer.nearColorTexture,
      this.renderer.midColorTexture,
      this.renderer.farColorTexture,
      this.renderer.nearDepthTexture,
      this.renderer.midDepthTexture,
      this.renderer.farDepthTexture,
      context.mainSceneTexture,
    );

    PostFXPass.clearTexture(this.renderer.core.device, this.renderer.orbitsTexture);
    this.renderer.postfxPass.run(encoder, context, theme, frameData.rawState.bodies ?? []);
  }
}



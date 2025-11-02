import type { Theme } from '@shared/types';
import type { RenderContext } from '../types';
import type { IRenderPipeline } from './base';
import { PostFXPass } from '../passes/postfxPass';
import { SceneRenderPass } from '../passes/sceneRenderPass';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { RenderPayload, ShipRelativePayload } from '@client/views/types';
import { CameraMode } from '@client/state';

export class ShipRelativePipeline implements IRenderPipeline {
  private sceneRenderPass: SceneRenderPass;
  private postfxPass: PostFXPass;

  constructor() {
    this.sceneRenderPass = new SceneRenderPass();
    this.postfxPass = new PostFXPass();
  }

  public async initialize(core: WebGPUCore, scene: Scene): Promise<void> {
    await this.sceneRenderPass.initialize(core, scene);
    await this.postfxPass.initialize(core, scene);
  }

  public onResize(_size: { width: number; height: number }, _core: WebGPUCore): void {
    // No intermediate textures to resize anymore.
  }

  render(encoder: GPUCommandEncoder, context: RenderContext, frameData: RenderPayload, theme: Theme): void {
    if (frameData.cameraMode !== CameraMode.SHIP_RELATIVE) return;

    if (frameData.dominantLight && frameData.worldCameraEye) {
      this.updateSceneUniforms(context.core, frameData, context);
      this.sceneRenderPass.run(encoder, context);
    } else {
      // Handle case with no light
      this.updateSceneUniforms(context.core, frameData, context);
      this.sceneRenderPass.run(encoder, context);
    }

    PostFXPass.clearTexture(context.core.device, context.orbitsTexture);
    this.postfxPass.run(encoder, context, theme, frameData.rawState.bodies ?? []);
  }

  private updateSceneUniforms(
    core: WebGPUCore,
    frameData: ShipRelativePayload,
    context: RenderContext
  ) {
    const light = frameData.dominantLight;
    const emissive = light?.emissive ?? [1, 1, 1];
    const intensity = 15.0;

    // Light direction is from world origin towards the dominant light source,
    // since all objects are now in a camera-relative frame.
    let lightDir: [number, number, number] = [
        light.position[0] - frameData.worldCameraEye[0],
        light.position[1] - frameData.worldCameraEye[1],
        light.position[2] - frameData.worldCameraEye[2],
    ];
    const lightDirLen = Math.hypot(lightDir[0], lightDir[1], lightDir[2]);
    if (lightDirLen > 0) {
      lightDir = [lightDir[0] / lightDirLen, lightDir[1] / lightDirLen, lightDir[2] / lightDirLen];
    }

    const len = Math.hypot(emissive[0], emissive[1], emissive[2]);
    const finalLightColor: [number, number, number] = len > 0
      ? [
          (emissive[0] / len) * intensity,
          (emissive[1] / len) * intensity,
          (emissive[2] / len) * intensity,
        ]
      : [0, 0, 0];
    const bufferData = new Float32Array(12);
    bufferData.set(lightDir, 0);
    bufferData[3] = 0.0;
    bufferData.set(finalLightColor, 4);
    bufferData[7] = frameData.debugTierView;
    bufferData[8] = 1.0; // tierScale is now always 1.0
    bufferData[9] = frameData.showAtmosphere ? 1.0 : 0.0;
    core.device.queue.writeBuffer(context.sceneUniformBuffer!, 0, bufferData);
  }

  public prepare(_frameData: RenderPayload, _renderer: any): void {
    // No-op for this pipeline
  }
}



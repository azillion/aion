import type { Camera } from '../camera';
import type { Scene } from './scene';
import type { WebGPUCore } from './core';

export interface RenderContext {
  core: WebGPUCore;
  scene: Scene;
  camera: Camera;
  systemScale: number;
  textureSize: { width: number, height: number };
  
  // For post-processing
  sourceTexture: GPUTexture;
  destinationTexture: GPUTexture;

  // Render targets
  mainSceneTexture: GPUTexture;
  orbitsTexture: GPUTexture;
  
  // Uniforms
  themeUniformBuffer: GPUBuffer;
  lastDeltaTime: number;
}

export interface IRenderPass {
  initialize(core: WebGPUCore, scene: Scene): void;
  run(encoder: GPUCommandEncoder, context: RenderContext, ...args: any[]): void;
  onResize?(size: { width: number, height: number }, core: WebGPUCore): void;
}

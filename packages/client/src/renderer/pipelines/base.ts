import type { FrameData, Theme } from '@shared/types';
import type { RenderContext } from '../types';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';

export interface IRenderPipeline {
  initialize(core: WebGPUCore, scene: Scene): Promise<void>;
  onResize(size: { width: number; height: number }, core: WebGPUCore): void;
  prepare(frameData: FrameData, renderer: any): void;
  render(encoder: GPUCommandEncoder, context: RenderContext, frameData: FrameData, theme: Theme): void;
}



import type { Theme } from '@shared/types';
import type { RenderContext } from '../types';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { RenderPayload } from '@client/orchestration/types';

export interface IRenderPipeline {
  initialize(core: WebGPUCore, scene: Scene): Promise<void>;
  onResize(size: { width: number; height: number }, core: WebGPUCore): void;
  prepare(frameData: RenderPayload, renderer: any): void;
  render(encoder: GPUCommandEncoder, context: RenderContext, frameData: RenderPayload, theme: Theme): void;
}



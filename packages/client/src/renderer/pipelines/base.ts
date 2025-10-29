import type { FrameData, Theme } from '@shared/types';
import type { RenderContext } from '../types';

export interface IRenderPipeline {
  render(encoder: GPUCommandEncoder, context: RenderContext, frameData: FrameData, theme: Theme): void;
}



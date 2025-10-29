import type { FrameData, Theme } from '../../shared/types';
import type { RenderContext } from '../types';
import type { GlyphsPass } from '../passes/glyphsPass';
import type { MapPass } from '../passes/mapPass';
import type { OrbitsPass } from '../passes/orbitsPass';
import type { SOIPass } from '../passes/soiPass';
import type { IRenderPipeline } from './base';
import type { AppState } from '../../state';

export class SystemMapPipeline implements IRenderPipeline {
  constructor(
    private orbitsPass: OrbitsPass,
    private mapPass: MapPass,
    private soiPass: SOIPass,
    private glyphsPass: GlyphsPass,
    private appState: AppState
  ) {}

  render(encoder: GPUCommandEncoder, context: RenderContext, _frameData: FrameData, theme: Theme): void {
    const showOrbits = this.appState.showOrbits;
    if (showOrbits) {
      this.orbitsPass.run(encoder, context, theme);
    }
    this.mapPass.run(encoder, context, theme, showOrbits);
    if (showOrbits) {
      this.soiPass.run(encoder, context);
      this.glyphsPass.run(encoder, context);
    }
  }
}



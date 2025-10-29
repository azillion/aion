import type { FrameData, Theme } from '@shared/types';
import type { RenderContext } from '../types';
import { GlyphsPass } from '../passes/glyphsPass';
import { MapPass } from '../passes/mapPass';
import { OrbitsPass } from '../passes/orbitsPass';
import { SOIPass } from '../passes/soiPass';
import type { IRenderPipeline } from './base';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';

export class SystemMapPipeline implements IRenderPipeline {
  private orbitsPass: OrbitsPass;
  private mapPass: MapPass;
  private soiPass: SOIPass;
  private glyphsPass: GlyphsPass;

  constructor() {
    this.orbitsPass = new OrbitsPass();
    this.mapPass = new MapPass();
    this.soiPass = new SOIPass();
    this.glyphsPass = new GlyphsPass();
  }

  public async initialize(core: WebGPUCore, scene: Scene): Promise<void> {
    await this.orbitsPass.initialize(core, scene);
    await this.mapPass.initialize(core, scene);
    await this.soiPass.initialize(core, scene);
    await this.glyphsPass.initialize(core, scene);
  }

  public onResize(_size: { width: number; height: number }, _core: WebGPUCore): void {
    // No-op
  }

  public prepare(frameData: FrameData, renderer: any): void {
    const { showOrbits, unscaledBodiesForMap, systemScale } = frameData;
    if (!showOrbits) return;
    if (!unscaledBodiesForMap) return;
    const glyphData = this.orbitsPass.update(unscaledBodiesForMap, renderer.getScene(), renderer.getCore(), systemScale);
    this.glyphsPass.update(
      glyphData.periapsisPoints,
      glyphData.apoapsisPoints,
      glyphData.ascendingNodePoints,
      glyphData.descendingNodePoints
    );
    this.soiPass.update(unscaledBodiesForMap, renderer.getScene(), systemScale);
  }

  public clearOrbitHistory(): void {
    this.orbitsPass.clearAll();
  }

  render(encoder: GPUCommandEncoder, context: RenderContext, frameData: FrameData, theme: Theme): void {
    const { showOrbits } = frameData;
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



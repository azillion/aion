import type { Body, SystemState } from '@shared/types';
import { ReferenceFrame } from '../state';
import type { Camera } from '../camera';
import { Scene } from './scene';

export interface SystemMapData {
  bodiesToRender: Body[];
  renderScale: number;
  // This contains the unscaled, but reference-frame-transformed bodies
  unscaledBodiesForMap: Body[];
}

export class SystemMapManager {
  public prepare(
    systemState: SystemState,
    referenceFrame: ReferenceFrame,
    camera: Camera
  ): SystemMapData {
    let bodiesForMap = systemState.bodies;
    // --- Reference Frame Transformation ---
    const focusBodyForFrame = systemState.bodies.find(b => b.id === camera.focusBodyId);
    if (referenceFrame === ReferenceFrame.FOCUSED_BODY && focusBodyForFrame) {
      bodiesForMap = systemState.bodies.map(body => ({
        ...body,
        position: [
          body.position[0] - focusBodyForFrame.position[0],
          body.position[1] - focusBodyForFrame.position[1],
          body.position[2] - focusBodyForFrame.position[2],
        ] as [number, number, number],
      }));
    }

    const finalRenderScale = Scene.calculateRenderScale(bodiesForMap, camera.focusBodyId);
    const bodiesToRender = bodiesForMap.map(b => ({
      ...b,
      position: [
        b.position[0] * finalRenderScale,
        b.position[1] * finalRenderScale,
        b.position[2] * finalRenderScale,
      ] as [number, number, number],
    }));

    return {
      bodiesToRender,
      renderScale: finalRenderScale,
      unscaledBodiesForMap: bodiesForMap,
    };
  }
}



import type { IRendererAccessor, IViewManager, RenderPayload, SystemMapPayload } from './types';
import type { SystemState } from '@shared/types';
import type { InputManager } from '@client/input';
import { SystemMapManager } from '@client/renderer/systemMapManager';
import { CameraMode } from '@client/state';
import { projectWorldToScreen } from '@client/renderer/projection';

export class SystemMapViewManager implements IViewManager {
  private systemMapManager: SystemMapManager;

  constructor() {
    this.systemMapManager = new SystemMapManager();
  }

  public prepare(
    systemState: SystemState,
    cameraManager: any,
    input: InputManager,
    viewport: { width: number, height: number },
    renderer: IRendererAccessor,
  ): RenderPayload {
    const camera = cameraManager.getCamera();
    const mapData = this.systemMapManager.prepare(systemState, renderer.state.referenceFrame, camera);
    const bodiesToRender = mapData.bodiesToRender;
    const renderScale = mapData.renderScale;

    cameraManager.update(CameraMode.SYSTEM_MAP, { bodies: mapData.unscaledBodiesForMap, scale: mapData.renderScale, viewport: viewport, vfov: 25.0, referenceFrame: renderer.state.referenceFrame });

    const unscaledBodiesForMap = mapData.unscaledBodiesForMap;

    // Handle click selection
    if (input.clicked) {
      let closestBodyId: string | null = null;
      let closestDist2 = Number.POSITIVE_INFINITY;
      for (const b of bodiesToRender) {
        const screenPos = projectWorldToScreen(b.position as [number,number,number], camera, viewport);
        if (!screenPos) continue;
        const dx = screenPos.x - input.mouseX;
        const dy = screenPos.y - input.mouseY;
        const d2 = dx * dx + dy * dy;
        if (d2 < closestDist2) {
          closestDist2 = d2;
          closestBodyId = b.id;
        }
      }
      if (closestBodyId) {
        camera.focusBodyId = closestBodyId;
        camera.pendingFrame = true;
        if (renderer.ui) {
          const bodyName = systemState.bodies.find(b => b.id === closestBodyId)?.name;
          if (bodyName) {
            renderer.ui.setFocus(bodyName);
          }
        }
      }
    }

    const payload: SystemMapPayload = {
      rawState: systemState,
      bodiesToRender,
      camera,
      systemScale: renderScale,
      viewport,
      deltaTime: 0, // Will be filled in by app loop
      cameraMode: CameraMode.SYSTEM_MAP,
      playerShipId: renderer.state.playerShipId,
      showOrbits: renderer.state.showOrbits,
      unscaledBodiesForMap,
    };
    return payload;
  }
}



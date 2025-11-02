import type { Body, SystemState } from '@shared/types';
import type { Camera } from '@client/camera/camera';
import type { InputManager } from '@client/input';
import { CameraMode } from '@client/state';
import type { Scene } from '@client/renderer/scene';
import type { UI } from '@client/ui';

// An accessor interface to break dependency cycles between ViewManagers and the full Renderer.
export interface IRendererAccessor {
  state: {
    playerShipId: string | null;
    showOrbits: boolean;
    showAtmosphere: boolean;
    referenceFrame: any; // Keep `any` for now for simplicity
  };
  ui: UI | null;
  getScene(): Scene;
}

// Base payload with common properties
interface BasePayload {
  rawState: SystemState;
  camera: Camera;
  viewport: { width: number; height: number };
  deltaTime: number;
  playerShipId: string | null;
}

// Payload specific to the ShipRelative view
export interface ShipRelativePayload extends BasePayload {
  cameraMode: CameraMode.SHIP_RELATIVE;
  worldCameraEye: [number, number, number];
  dominantLight: Body;
  showAtmosphere: boolean;
  // This list is camera-relative and prepared specifically for the HUD
  bodiesToRender: Body[];
}

export type RenderPayload = ShipRelativePayload; // Only one payload type for now

export interface IViewManager {
  /**
   * Takes the raw simulation state and prepares a tailored RenderPayload
   * for the current view mode. This includes camera updates, culling,
   * scaling, and any other view-specific logic.
   */
  prepare(
    state: SystemState,
    cameraManager: any, // Using 'any' for now to ease refactoring
    input: InputManager,
    viewport: { width: number, height: number },
    renderer: IRendererAccessor,
  ): RenderPayload;
}



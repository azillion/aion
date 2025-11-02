import type { IRendererAccessor, IViewManager, RenderPayload, ShipRelativePayload } from './types';
import type { SystemState, Vec3, Body } from '@shared/types';
import type { InputManager } from '@client/input';
import { CameraMode } from '@client/state';
import type { Scene } from '@client/renderer/scene';
import { SceneDataProcessor } from './sceneDataProcessor';

export class ShipRelativeViewManager implements IViewManager {
  private sceneDataProcessor: SceneDataProcessor;

  constructor() {
    this.sceneDataProcessor = new SceneDataProcessor();
  }

  public prepare(
    systemState: SystemState,
    cameraManager: any,
    input: InputManager,
    viewport: { width: number, height: number },
    renderer: IRendererAccessor,
  ): RenderPayload {
    const scene: Scene = renderer.getScene();
    const camera = cameraManager.getCamera();
    const state = renderer.state;

    const playerShip = systemState.ships.find(s => s.body.id === state.playerShipId);

    // 1) Update world-space camera from controller FIRST.
    cameraManager.update(CameraMode.SHIP_RELATIVE, { playerShip, keys: input.keys, viewport });

    const worldCameraEye = [...camera.eye] as Vec3; // Capture true f64 (JS number) camera position.
    
    // 2) Prepare a single list of all bodies with their true f64 world positions.
    // No more tiering or scaling.
    const bodiesToRender: Body[] = systemState.bodies.filter(b => b.id !== state.playerShipId);
    const shadowCasters: Body[] = [];
    systemState.bodies.forEach(body => {
        if (body.mass > 1e22) {
            shadowCasters.push({ ...body }); // Pass absolute f64 world position
        }
    });

    // 3) Serialize the data for the GPU.
    const processedData = this.sceneDataProcessor.process(bodiesToRender);

    // 6) Move the main camera to the origin for this frame's rendering.
    camera.eye = [0, 0, 0];

    // 5) Compute look_at using camera-relative data for precision.
    const relativeBodiesForLookAt = bodiesToRender.map(b => ({
        ...b,
        position: [b.position[0] - worldCameraEye[0], b.position[1] - worldCameraEye[1], b.position[2] - worldCameraEye[2]] as Vec3
    }));
    cameraManager.updateLookAt(camera, { keys: input.keys, relativeBodies: relativeBodiesForLookAt, playerShip });
    
    // 6) Upload data to GPU buffers.
    scene.updateSceneObjects(processedData.data, processedData.count);
    // Serialize shadow casters using the same f64-split logic.
    const processedShadowData = this.sceneDataProcessor.process(shadowCasters);
    scene.updateShadowCasters(processedShadowData.data, processedShadowData.count);

    // Determine dominant light by emissive energy
    let dominantLight = systemState.bodies[0];
    let maxEmissiveEnergy = 0.0;
    systemState.bodies.forEach(body => {
      if (body.emissive) {
        const energy = body.emissive[0] + body.emissive[1] + body.emissive[2];
        if (energy > maxEmissiveEnergy) {
          maxEmissiveEnergy = energy;
          dominantLight = body;
        }
      }
    });
    (camera as any).dominantLight = dominantLight;

    // Handle targeting key
    if (input.wasPressed('KeyN')) {
      if (playerShip) {
        let closestBodyId: string | null = null;
        let closestDistSq = Infinity;
        for (const body of systemState.bodies) {
          if (body.id === state.playerShipId) continue;
          const dx = body.position[0] - playerShip.body.position[0];
          const dy = body.position[1] - playerShip.body.position[1];
          const dz = body.position[2] - playerShip.body.position[2];
          const distSq = dx*dx + dy*dy + dz*dz;
          if (distSq < closestDistSq) {
            closestDistSq = distSq;
            closestBodyId = body.id;
          }
        }
        if (closestBodyId) {
          camera.focusBodyId = closestBodyId;
          camera.pendingFrame = true;
        }
      }
    }

    const payload: ShipRelativePayload = {
      rawState: systemState,
      camera,
      viewport,
      deltaTime: 0, // Will be filled in by app loop
      cameraMode: CameraMode.SHIP_RELATIVE,
      playerShipId: state.playerShipId,
      dominantLight: dominantLight,
      worldCameraEye: worldCameraEye,
      debugTierView: state.debugTierView,
      showAtmosphere: state.showAtmosphere,
      // Pass the specially-prepared HUD body list
      bodiesToRender: bodiesToRender, // Pass world-space bodies to HUD
    };
    return payload;
  }
}



import type { IRendererAccessor, IViewManager, RenderPayload, ShipRelativePayload } from './types';
import type { SystemState, Vec3 } from '@shared/types';
import type { InputManager } from '@client/input';
import { TierManager } from '@client/renderer/tierManager';
import { CameraMode } from '@client/state';
import type { Scene } from '@client/renderer/scene';
import { SceneDataProcessor } from './sceneDataProcessor';

export class ShipRelativeViewManager implements IViewManager {
  private tierManager: TierManager;
  private sceneDataProcessor: SceneDataProcessor;

  constructor() {
    this.tierManager = new TierManager();
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
    const renderScale = 1.0;

    const playerShip = systemState.ships.find(s => s.body.id === state.playerShipId);

    // 1) Update world-space camera from controller FIRST.
    cameraManager.update(CameraMode.SHIP_RELATIVE, { playerShip, keys: input.keys, viewport });

    const shipCameraEyeWorld = [...camera.eye] as Vec3; // Capture true f64 (JS number) camera position.
    
    // Tiering and serialization (CPU, f64 -> f32-safe, zero-padded per tier)
    const tieredScene = this.tierManager.build(systemState, shipCameraEyeWorld, state.playerShipId);
    const processedData = this.sceneDataProcessor.process(tieredScene);

    // 4) HUD/selection should use unscaled, camera-relative positions.
    const bodiesForTiers = systemState.bodies.filter(b => b.id !== state.playerShipId);
    const bodiesForHUD = bodiesForTiers
      .map(b => {
        const HUD_RENDER_DISTANCE = 5000; // keep within stable f32 range for projection
        const p: [number, number, number] = [
          b.position[0] - shipCameraEyeWorld[0],
          b.position[1] - shipCameraEyeWorld[1],
          b.position[2] - shipCameraEyeWorld[2],
        ];
        const dist = Math.hypot(p[0], p[1], p[2]);
        if (dist > HUD_RENDER_DISTANCE) {
          const inv = 1.0 / dist;
          p[0] = p[0] * inv * HUD_RENDER_DISTANCE;
          p[1] = p[1] * inv * HUD_RENDER_DISTANCE;
          p[2] = p[2] * inv * HUD_RENDER_DISTANCE;
        }
        return { ...b, position: p };
      });

    // 6) Move the main camera to the origin for this frame's rendering.
    camera.eye = [0, 0, 0];

    // Tiered data already computed above

    // 7) Compute look_at using stable, camera-relative data.
    cameraManager.updateLookAt(camera, { keys: input.keys, relativeBodies: bodiesForHUD, playerShip });
    
    // Upload tier data to GPU buffers for ship-relative rendering
    const sceneTyped: Scene = scene as unknown as Scene;
    sceneTyped.updateTiers(
      processedData.nearData,
      processedData.midData,
      processedData.farData,
      processedData.nearCount,
      processedData.midCount,
      processedData.farCount
    );
    // Upload global shadow casters (unscaled camera-relative positions)
    sceneTyped.updateShadowCasters(tieredScene.shadowCasters);

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
      worldCameraEye: shipCameraEyeWorld,
      debugTierView: state.debugTierView,
      showAtmosphere: state.showAtmosphere,
      // Pass the specially-prepared HUD body list
      bodiesToRender: bodiesForHUD,
    };
    return payload;
  }
}



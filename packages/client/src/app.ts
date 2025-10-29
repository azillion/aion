import type { Authority } from './authority';
import type { InputManager } from './input';
import type { Renderer } from './renderer/index';
import type { AppState } from './state';
import type { UI } from './ui';
import type { HUDManager } from './hud';
import { Scene } from './renderer/scene';
import { CameraManager } from './camera/manager';
import type { Body, FrameData, Ship } from '@shared/types';
import { PLANETARY_SOI_RADIUS_MULTIPLIER } from '@shared/constants';
import { TierManager } from './renderer/tierManager';
import { SystemMapManager } from './renderer/systemMapManager';
import { CameraMode } from './state';
import { projectWorldToScreen } from './renderer/projection';
import { type OrbitGlyphData } from './renderer/passes/orbitsPass';
import Stats from 'stats.js';

export class App {
  private readonly authority: Authority;
  private readonly renderer: Renderer;
  private readonly state: AppState;
  private readonly input: InputManager;
  private readonly ui: UI;
  private readonly hud: HUDManager;
  private scene!: Scene;
  private cameraManager: CameraManager;
  private tierManager: TierManager;
  private systemMapManager: SystemMapManager;
  private lastFrameTime: number = 0;
  private fpsAccumulator = 0;
  private fpsFrameCount = 0;
  private fpsLastReport = 0;
  private stats: Stats | null = null;

  

  constructor(
    renderer: Renderer,
    authority: Authority,
    state: AppState,
    input: InputManager,
    ui: UI,
    hud: HUDManager,
    cameraManager: CameraManager
  ) {
    this.renderer = renderer;
    this.authority = authority;
    this.state = state;
    this.input = input;
    this.ui = ui;
    this.hud = hud;
    this.cameraManager = cameraManager;
    this.tierManager = new TierManager();
    this.systemMapManager = new SystemMapManager();
  }

  public async start(): Promise<void> {
    await this.renderer.initialize(this.authority);
    this.scene = this.renderer.getScene();
    // Initialize Stats.js overlay (FPS panel)
    this.stats = new Stats();
    this.stats.showPanel(0);
    document.body.appendChild(this.stats.dom);
    requestAnimationFrame(this.updateLoop);
  }

  private updateLoop = async (time: number): Promise<void> => {
    if (this.stats) this.stats.begin();
    if (this.lastFrameTime === 0) this.lastFrameTime = time;
    const deltaTime = (time - this.lastFrameTime) / 1000.0;
    this.lastFrameTime = time;

    // Accumulate FPS and update UI ~2 times per second
    this.fpsAccumulator += deltaTime;
    this.fpsFrameCount += 1;
    if ((time - this.fpsLastReport) > 500) {
      const avgDelta = this.fpsAccumulator / Math.max(1, this.fpsFrameCount);
      const fps = avgDelta > 0 ? (1 / avgDelta) : 0;
      if (this.ui && typeof (this.ui as any).setFps === 'function') {
        (this.ui as any).setFps(fps);
      }
      this.fpsAccumulator = 0;
      this.fpsFrameCount = 0;
      this.fpsLastReport = time;
    }

    const clampedDeltaTime = Math.min(deltaTime, 1 / 20.0);
    await this.authority.tick(clampedDeltaTime, { deltaX: this.input.deltaX, deltaY: this.input.deltaY, keys: this.input.keys });
    
    const systemState = await this.authority.query();

    // --- Reference Frame Management ---
    const worldCameraEye = [...this.cameraManager.getCamera().eye];
    let shipWorldCameraEye: [number, number, number] | null = null;
    const currentReferenceBody = systemState.bodies.find(b => b.id === this.state.referenceBodyId) ?? systemState.bodies.find(b => b.id === 'sol');
    let newReferenceBody = currentReferenceBody;
    let closestDistSq = Infinity;

    systemState.bodies.forEach(body => {
      if (body.mass < 1e22) return; // Only consider sufficiently massive bodies
      const dx = body.position[0] - worldCameraEye[0];
      const dy = body.position[1] - worldCameraEye[1];
      const dz = body.position[2] - worldCameraEye[2];
      const distSq = dx*dx + dy*dy + dz*dz;
      const soiRadius = body.radius * PLANETARY_SOI_RADIUS_MULTIPLIER;
      if (distSq < (soiRadius * soiRadius) && distSq < closestDistSq) {
        closestDistSq = distSq;
        newReferenceBody = body;
      }
    });
    if (closestDistSq === Infinity && this.state.referenceBodyId !== 'sol') {
      newReferenceBody = systemState.bodies.find(b => b.id === 'sol');
    }
    if (newReferenceBody && this.state.referenceBodyId !== newReferenceBody.id) {
      console.log(`Switching reference frame to: ${newReferenceBody.name}`);
      this.state.referenceBodyId = newReferenceBody.id;
    }
    // --- End Reference Frame Management ---

    let bodiesToRender: Body[] = systemState.bodies;
    let renderScale: number;
    let camera = this.cameraManager.getCamera();

    if (this.state.cameraMode === CameraMode.SYSTEM_MAP) {
      const mapData = this.systemMapManager.prepare(systemState, this.state.referenceFrame, camera);
      bodiesToRender = mapData.bodiesToRender;
      renderScale = mapData.renderScale;
      this.cameraManager.update(this.state.cameraMode, { bodies: mapData.unscaledBodiesForMap, scale: mapData.renderScale, viewport: this.renderer.getTextureSize(), vfov: 25.0, referenceFrame: this.state.referenceFrame });
      if (this.state.showOrbits) {
        // Recompute orbit geometry at correct scale
        const glyphData: OrbitGlyphData = this.renderer.getOrbitsPass().update(mapData.unscaledBodiesForMap, this.scene, this.renderer.getCore(), mapData.renderScale);
        this.renderer.getGlyphsPass().update(
          glyphData.periapsisPoints,
          glyphData.apoapsisPoints,
          glyphData.ascendingNodePoints,
          glyphData.descendingNodePoints
        );
        this.renderer.getSoiPass().update(mapData.unscaledBodiesForMap, this.scene, mapData.renderScale);
      }
    } else if (this.state.cameraMode === CameraMode.SHIP_RELATIVE) {
      renderScale = 1.0;
      const playerShip = systemState.bodies.find(b => b.id === this.state.playerShipId) as Ship | undefined;

      // 1) Update world-space camera from controller FIRST.
      this.cameraManager.update(this.state.cameraMode, { playerShip, keys: this.input.keys });
      const shipCameraEyeWorld = [...camera.eye]; // Capture true f64 (JS number) camera position.
      const tieredScene = this.tierManager.build(systemState, shipCameraEyeWorld as [number, number, number], this.state.playerShipId);

      // 4) HUD/selection should use unscaled, camera-relative positions.
      // Provide a separate list for HUD logic while renderer uses tier buffers directly.
      bodiesToRender = systemState.bodies
        .filter(b => b.id !== this.state.playerShipId)
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

      shipWorldCameraEye = shipCameraEyeWorld as [number, number, number];

      // 6) Move the main camera to the origin for this frame's rendering.
      camera.eye = [0, 0, 0];

      // 7) Compute look_at using stable, camera-relative data.
      // We pass the combined list here, as the controller might need to find a target.
      this.cameraManager.updateLookAt(camera, { keys: this.input.keys, relativeBodies: bodiesToRender, playerShip });
      // Upload tier data to GPU buffers for ship-relative rendering
      this.scene.updateTiers(tieredScene.near, tieredScene.mid, tieredScene.far);
      // Upload global shadow casters (unscaled camera-relative positions)
      this.scene.updateShadowCasters(tieredScene.shadowCasters);
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
    } else {
      renderScale = 1.0;
      this.cameraManager.update(this.state.cameraMode, {});
    }

    // Finalize camera and write shared camera buffer
    camera.updateViewMatrix();
    this.renderer.writeCameraBuffer(camera);

    const sceneScale = 1.0; // All scaling is now pre-applied in the tier sort.
    if (this.state.cameraMode === CameraMode.SYSTEM_MAP) {
      this.scene.updateMapBuffer(
        systemState,
        camera,
        bodiesToRender.filter(b => b.id !== this.state.playerShipId),
        sceneScale,
        true
      );
    }
    // Do not update orbits in ship view; rings are hidden there

    const textureSize = this.renderer.getTextureSize();

    const frameData: FrameData | null = textureSize ? {
      rawState: systemState,
      bodiesToRender,
      camera,
      systemScale: renderScale,
      viewport: textureSize,
      deltaTime,
      cameraMode: this.state.cameraMode,
      playerShipId: this.state.playerShipId,
      dominantLight: (camera as any).dominantLight,
      worldCameraEye: (shipWorldCameraEye ?? worldCameraEye) as [number, number, number],
      debugTierView: this.state.debugTierView,
    } : null;

    // --- Handle Ship-Relative Targeting Keys ---
    if (frameData && this.state.cameraMode === CameraMode.SHIP_RELATIVE) {
      // KeyN: Target nearest body (camera-relative)
      if (this.input.wasPressed('KeyN')) {
        const playerShip = systemState.bodies.find(b => b.id === this.state.playerShipId);
        if (playerShip) {
          let closestBodyId: string | null = null;
          let closestDistSq = Infinity;
          for (const body of systemState.bodies) {
            if (body.id === this.state.playerShipId) continue;
            const dx = body.position[0] - playerShip.position[0];
            const dy = body.position[1] - playerShip.position[1];
            const dz = body.position[2] - playerShip.position[2];
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
    }

    // Handle click selection in System Map using frameData
    if (frameData && this.state.cameraMode === CameraMode.SYSTEM_MAP && this.input.clicked) {
      const viewport = frameData.viewport;
      let closestBodyId: string | null = null;
      let closestDist2 = Number.POSITIVE_INFINITY;
      for (const b of frameData.bodiesToRender) {
        const screenPos = projectWorldToScreen(b.position as [number,number,number], camera, viewport);
        if (!screenPos) continue;

        const dx = screenPos.x - this.input.mouseX;
        const dy = screenPos.y - this.input.mouseY;
        const d2 = dx * dx + dy * dy;
        if (d2 < closestDist2) {
          closestDist2 = d2;
          closestBodyId = b.id;
        }
      }
      if (closestBodyId) {
        camera.focusBodyId = closestBodyId;
        camera.pendingFrame = true;
        if (this.ui) {
          const bodyName = systemState.bodies.find(b => b.id === closestBodyId)?.name;
          if (bodyName) {
            this.ui.setFocus(bodyName);
          }
        }
      }
    }

    this.renderer.render(frameData ?? { rawState: systemState, bodiesToRender, camera, systemScale: renderScale, viewport: { width: 0, height: 0 }, deltaTime, cameraMode: this.state.cameraMode, playerShipId: this.state.playerShipId });

    if (frameData && this.hud) {
      if (this.state.showHUD) {
        this.hud.draw(frameData);
      } else {
        this.hud.clear();
      }
    }

    this.input.tick();
    if (this.stats) this.stats.end();
    requestAnimationFrame(this.updateLoop);
  };
}



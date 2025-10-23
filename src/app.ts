import type { Authority } from './authority/authority';
import type { InputManager } from './input';
import type { Renderer } from './renderer/index';
import type { AppState } from './state';
import type { UI } from './ui';
import type { HUDManager } from './hud';
import { Scene } from './renderer/scene';
import { CameraManager } from './camera/manager';
import type { Body, FrameData, Ship } from './shared/types';
import { CameraMode, ReferenceFrame } from './state';
import { projectWorldToScreen } from './renderer/projection';
import { type OrbitGlyphData } from './renderer/passes/orbitsPass';

export class App {
  private readonly authority: Authority;
  private readonly renderer: Renderer;
  private readonly state: AppState;
  private readonly input: InputManager;
  private readonly ui: UI;
  private readonly hud: HUDManager;
  private scene!: Scene;
  private cameraManager: CameraManager;
  private lastFrameTime: number = 0;

  

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
  }

  public async start(): Promise<void> {
    await this.renderer.initialize(this.authority);
    this.scene = this.renderer.getScene();
    requestAnimationFrame(this.updateLoop);
  }

  private updateLoop = async (time: number): Promise<void> => {
    if (this.lastFrameTime === 0) this.lastFrameTime = time;
    const deltaTime = (time - this.lastFrameTime) / 1000.0;
    this.lastFrameTime = time;

    const clampedDeltaTime = Math.min(deltaTime, 1 / 20.0);
    await this.authority.tick(clampedDeltaTime, { deltaX: this.input.deltaX, deltaY: this.input.deltaY, keys: this.input.keys });
    
    const systemState = await this.authority.query();

    let bodiesToRender: Body[] = systemState.bodies;
    let renderScale: number;
    let camera = this.cameraManager.getCamera();

    if (this.state.cameraMode === CameraMode.SYSTEM_MAP) {
      let bodiesForMap = systemState.bodies;
      // --- Reference Frame Transformation ---
      const focusBodyForFrame = systemState.bodies.find(b => b.id === camera.focusBodyId);
      if (this.state.referenceFrame === ReferenceFrame.FOCUSED_BODY && focusBodyForFrame) {
        bodiesForMap = systemState.bodies.map(body => ({
          ...body,
          position: [
            body.position[0] - focusBodyForFrame.position[0],
            body.position[1] - focusBodyForFrame.position[1],
            body.position[2] - focusBodyForFrame.position[2],
          ] as [number, number, number],
        }));
      }

      renderScale = Scene.calculateRenderScale(bodiesForMap, camera.focusBodyId);
      this.cameraManager.update(this.state.cameraMode, { bodies: bodiesForMap, scale: renderScale, viewport: this.renderer.getTextureSize(), vfov: 25.0, referenceFrame: this.state.referenceFrame });
      bodiesToRender = bodiesForMap.map(b => ({
        ...b,
        position: [
          b.position[0] * renderScale,
          b.position[1] * renderScale,
          b.position[2] * renderScale,
        ] as [number, number, number],
      }));
      if (this.state.showOrbits) {
        // Recompute orbit geometry at correct scale
        const glyphData: OrbitGlyphData = this.renderer.getOrbitsPass().update(bodiesForMap, this.scene, this.renderer.getCore(), renderScale);
        this.renderer.getGlyphsPass().update(
          glyphData.periapsisPoints,
          glyphData.apoapsisPoints,
          glyphData.ascendingNodePoints,
          glyphData.descendingNodePoints
        );
        this.renderer.getSoiPass().update(bodiesForMap, this.scene, renderScale);
      }
    } else if (this.state.cameraMode === CameraMode.SHIP_RELATIVE) {
      renderScale = 1.0;
      const playerShip = systemState.bodies.find(b => b.id === this.state.playerShipId) as Ship | undefined;
      // 1) Update world-space camera from controller (no look_at here)
      this.cameraManager.update(this.state.cameraMode, { playerShip, keys: this.input.keys });

      // 2) Build camera-relative bodies using the (world) camera.eye
      bodiesToRender = systemState.bodies.map(b => ({
        ...b,
        position: [
          b.position[0] - camera.eye[0],
          b.position[1] - camera.eye[1],
          b.position[2] - camera.eye[2],
        ] as [number, number, number],
      }));

      // 3) Move the main camera to origin for the rest of this frame
      camera.eye = [0, 0, 0];

      // 4) Compute look_at using stable, camera-relative data
      this.cameraManager.updateLookAt(camera, { keys: this.input.keys, relativeBodies: bodiesToRender, playerShip });
    } else {
      renderScale = 1.0;
      this.cameraManager.update(this.state.cameraMode, {});
    }

    // Finalize camera and write shared camera buffer
    camera.updateViewMatrix();
    this.renderer.writeCameraBuffer(camera);

    const sceneScale = this.state.cameraMode === CameraMode.SYSTEM_MAP ? 1.0 : renderScale;
    this.scene.update(
      systemState,
      camera,
      bodiesToRender.filter(b => b.id !== this.state.playerShipId),
      sceneScale,
      true
    );
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
    } : null;

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
      this.hud.draw(frameData);
    }

    this.input.tick();
    requestAnimationFrame(this.updateLoop);
  };
}



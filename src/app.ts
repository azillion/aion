import type { Authority } from './authority/authority';
import type { InputManager } from './input';
import type { Renderer } from './renderer/index';
import type { AppState } from './state';
import type { UI } from './ui';
import type { HUDManager } from './hud';
import { Scene } from './renderer/scene';
import { CameraManager } from './camera/manager';
import type { Body, Ship, SystemState, Vec3 } from './shared/types';
import { CameraMode, ReferenceFrame } from './state';
import { vec4 } from 'gl-matrix';

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
  private lastSystemState?: SystemState;

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
    this.lastSystemState = systemState;
    
    if (systemState.bodies.length !== this.scene.lastKnownBodyCount) {
        this.scene.recreateSpheresBuffer(systemState.bodies.length);
        this.renderer.getComputePass().recreatePipeline(systemState.bodies.length);
    }

    let bodiesToRender: Body[] = systemState.bodies;
    let renderScale: number;
    const camera = this.cameraManager.getCamera();

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
      bodiesToRender = bodiesForMap;
      if (this.state.showOrbits) {
        // Recompute orbit geometry at correct scale
        this.renderer.getOrbitsPass().update(bodiesForMap, this.scene, this.renderer.getCore(), renderScale);
        this.renderer.getGlyphsPass().update(
          this.renderer.getOrbitsPass().periapsisPoints,
          this.renderer.getOrbitsPass().apoapsisPoints,
          this.renderer.getOrbitsPass().ascendingNodePoints,
          this.renderer.getOrbitsPass().descendingNodePoints
        );
        this.renderer.getSoiPass().update(bodiesForMap, this.scene, renderScale);
      }
    } else if (this.state.cameraMode === CameraMode.SHIP_RELATIVE) {
      renderScale = 1.0;
      const playerShip = systemState.bodies.find(b => b.id === this.state.playerShipId) as Ship | undefined;
      const earth = systemState.bodies.find(b => b.id.toLowerCase() === 'earth');
      this.cameraManager.update(this.state.cameraMode, { playerShip, targetBody: earth });
      bodiesToRender = systemState.bodies.filter(b => b.id !== this.state.playerShipId);
    } else {
      renderScale = 1.0;
      this.cameraManager.update(this.state.cameraMode, {});
    }

    // Finalize camera and write shared camera buffer
    camera.updateViewMatrix();
    this.renderer.writeCameraBuffer(camera);

    this.scene.update(
      systemState,
      camera,
      bodiesToRender,
      renderScale,
      this.state.cameraMode === CameraMode.SYSTEM_MAP
    );
    // Do not update orbits in ship view; rings are hidden there

    const textureSize = this.renderer.getTextureSize();
    // Handle click selection in System Map
    if (this.state.cameraMode === CameraMode.SYSTEM_MAP && this.input.clicked && textureSize) {
      const viewport = textureSize;
      const viewProj = camera.viewProjectionMatrix as unknown as number[];
      let closestBodyId: string | null = null;
      let closestDist2 = Number.POSITIVE_INFINITY;
      // Recompute the transformed bodies for hit testing, same as render path
      const fb = systemState.bodies.find(b => b.id === camera.focusBodyId);
      const transformedForHit = (this.state.referenceFrame === ReferenceFrame.FOCUSED_BODY && fb)
        ? systemState.bodies.map(body => ({
            id: body.id,
            name: body.name,
            position: [
              body.position[0] - fb.position[0],
              body.position[1] - fb.position[1],
              body.position[2] - fb.position[2],
            ] as [number, number, number],
          }))
        : systemState.bodies;
      for (const b of transformedForHit) {
        const worldPos = vec4.fromValues(b.position[0] * renderScale, b.position[1] * renderScale, b.position[2] * renderScale, 1.0);
        const clip = vec4.transformMat4(vec4.create(), worldPos, viewProj as unknown as number[]);
        const w = clip[3];
        if (!Number.isFinite(w) || w === 0) continue;
        const ndcX = clip[0] / w;
        const ndcY = clip[1] / w;
        const px = (ndcX * 0.5 + 0.5) * viewport.width;
        const py = (-ndcY * 0.5 + 0.5) * viewport.height;
        const dx = px - this.input.mouseX;
        const dy = py - this.input.mouseY;
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
    
    this.renderer.render(camera, renderScale, deltaTime, this.lastSystemState);
    
    if (this.hud && textureSize) {
        let hudBodies = systemState.bodies;
        if (this.state.cameraMode === CameraMode.SYSTEM_MAP) {
        const fb2 = systemState.bodies.find(b => b.id === camera.focusBodyId);
          const bodiesForHud = (this.state.referenceFrame === ReferenceFrame.FOCUSED_BODY && fb2)
            ? systemState.bodies.map(b => ({
                ...b,
                position: [
                  b.position[0] - fb2.position[0],
                  b.position[1] - fb2.position[1],
                  b.position[2] - fb2.position[2],
                ] as [number, number, number],
              }))
            : systemState.bodies;
          const scaledBodies = bodiesForHud.map(b => ({ ...b, position: [ b.position[0] * renderScale, b.position[1] * renderScale, b.position[2] * renderScale ] as Vec3 }));
          hudBodies = scaledBodies;
        }
        this.hud.draw(hudBodies, camera, textureSize, this.state.cameraMode, this.state.playerShipId);
    }

    this.input.tick();
    requestAnimationFrame(this.updateLoop);
  };
}



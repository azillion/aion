import type { Authority } from './authority';
import type { InputManager } from './input';
import type { Renderer } from './renderer/index';
import type { AppState } from './state';
import type { UI } from './ui';
import type { HUDManager } from './hud';
import { CameraManager } from './camera/manager';
import { ShipRelativeViewManager } from './views/shipRelativeViewManager';
 
import Stats from 'stats.js';

export class App {
  private readonly authority: Authority;
  private readonly renderer: Renderer;
  private readonly state: AppState;
  private readonly input: InputManager;
  private readonly ui: UI;
  private readonly hud: HUDManager;
  private cameraManager: CameraManager;
  private shipRelativeViewManager: ShipRelativeViewManager;
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
    this.shipRelativeViewManager = new ShipRelativeViewManager();
  }

  public async start(): Promise<void> {
    await this.renderer.initialize(this.authority);
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
    await this.authority.tick(clampedDeltaTime, { deltaX: this.input.deltaX, deltaY: this.input.deltaY, keys: this.input.keys } as any);

    
    const systemState = await this.authority.query();
    // --- Prepare scene ---
    const textureSize = this.renderer.getTextureSize();
    if (!textureSize) {
      requestAnimationFrame(this.updateLoop);
      return;
    }

    const renderPayload = this.shipRelativeViewManager.prepare(systemState, this.cameraManager, this.input, textureSize, this.renderer as any);
    renderPayload.deltaTime = deltaTime; // Inject the final delta time

    // Finalize camera and write shared camera buffer
    const camera = this.cameraManager.getCamera();
    camera.updateViewMatrix();
    this.renderer.writeCameraBuffer(camera, renderPayload.worldCameraEye);

    this.renderer.prepare(renderPayload);
    this.renderer.render(renderPayload);

    if (this.hud) {
      if (this.state.showHUD) {
        this.hud.draw(renderPayload);
      } else {
        this.hud.clear();
      }
    }

    this.input.tick();
    if (this.stats) this.stats.end();
    requestAnimationFrame(this.updateLoop);
  };
}



import type { Authority } from '../authority/authority';
import type { InputManager } from '../input';
import { AppState, CameraMode, ReferenceFrame } from '../state';
import { CameraManager } from '../camera/manager';
import type { Camera } from '../camera';
import { HUDManager } from '../hud';
import { Galaxy } from '../galaxy';
import { WebGPUCore } from './core';
import { Scene } from './scene';
import type { RenderContext } from './types';
import { vec4, vec3 } from 'gl-matrix';
import type { UI } from '../ui';
import { ComputePass } from './passes/computePass';
import { GalaxyPass } from './passes/galaxyPass';
import { OrbitsPass } from './passes/orbitsPass';
import { PostFXPass } from './passes/postfxPass';
import { MapPass } from './passes/mapPass';
import { GlyphsPass } from './passes/glyphsPass';
import { SOIPass } from './passes/soiPass';
import { themes } from '../theme';
import { spectralResponses } from '../spectral';
import type { Body, Ship, Theme, Vec3 } from '../shared/types';

export class Renderer {
  public readonly authority: Authority;
  private readonly canvas: HTMLCanvasElement;
  private readonly state: AppState;
  private readonly input: InputManager;
  public ui: UI | null = null;
  
  private core!: WebGPUCore;
  private scene!: Scene;
  private cameraManager: CameraManager;
  private hud!: HUDManager;
  private galaxy: Galaxy;
  
  private computePass!: ComputePass;
  private galaxyPass!: GalaxyPass;
  private orbitsPass!: OrbitsPass;
  private postfxPass!: PostFXPass;
  private mapPass!: MapPass;
  private glyphsPass!: GlyphsPass;
  private soiPass!: SOIPass;

  private textureSize!: { width: number, height: number };
  private mainSceneTexture!: GPUTexture;
  private postFxTextureA!: GPUTexture;
  private postFxTextureB!: GPUTexture;
  private orbitsTexture!: GPUTexture;
  
  private themeUniformBuffer!: GPUBuffer;
  private currentThemeName: string = 'white';
  private currentResponseName: string = 'Visible (Y)';
  private lastDeltaTime: number = 1 / 60;
  private lastFrameTime: number = 0;
  private frameCount = 0;
  
  private lastSystemState?: import('../shared/types').SystemState;

  constructor(canvas: HTMLCanvasElement, authority: Authority, state: AppState, input: InputManager) {
    this.canvas = canvas;
    this.authority = authority;
    this.state = state;
    this.input = input;
    this.cameraManager = new CameraManager(this.canvas);
    this.galaxy = new Galaxy(1337, 20000, 200);

    const hudCanvas = document.getElementById('hud-canvas') as HTMLCanvasElement;
    if (hudCanvas) this.hud = new HUDManager(hudCanvas);
  }

  public getCamera(): Camera {
    return this.cameraManager.getCamera();
  }

  public async start() {
    this.core = await WebGPUCore.create(this.canvas);
    this.scene = new Scene(this.core);

    const initialSystemState = await this.authority.query();
    this.lastSystemState = initialSystemState;
    this.scene.initialize(initialSystemState, this.galaxy.stars);
    
    this.themeUniformBuffer = this.core.device.createBuffer({
      size: 80,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    this.computePass = new ComputePass();
    this.galaxyPass = new GalaxyPass();
    this.orbitsPass = new OrbitsPass();
    this.postfxPass = new PostFXPass();
    this.mapPass = new MapPass();
    this.glyphsPass = new GlyphsPass();
    this.soiPass = new SOIPass();
    
    this.computePass.initialize(this.core, this.scene);
    this.galaxyPass.initialize(this.core, this.scene);
    this.orbitsPass.initialize(this.core, this.scene);
    this.postfxPass.initialize(this.core, this.scene);
    this.mapPass.initialize(this.core, this.scene);
    this.glyphsPass.initialize(this.core, this.scene);
    this.soiPass.initialize(this.core, this.scene);

    this.handleResize();
    
    new ResizeObserver(this.handleResize).observe(this.canvas);
    requestAnimationFrame(this.updateLoop);
  }

  private handleResize = () => {
    const width = this.canvas.clientWidth;
    const height = this.canvas.clientHeight;
    if (width === 0 || height === 0) return;
    if (width === this.textureSize?.width && height === this.textureSize?.height) return;

    this.canvas.width = width;
    this.canvas.height = height;
    if (this.hud) {
        this.hud.getCanvas().width = width;
        this.hud.getCanvas().height = height;
    }
    this.textureSize = { width, height };
    this.core.context.configure({ device: this.core.device, format: this.core.presentationFormat });

    [this.mainSceneTexture, this.postFxTextureA, this.postFxTextureB, this.orbitsTexture]
        .forEach(tex => tex?.destroy());

    this.mainSceneTexture = this.core.device.createTexture({ size: this.textureSize, format: 'rgba16float', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING });
    this.postFxTextureA = this.core.device.createTexture({ size: this.textureSize, format: this.core.presentationFormat, usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC });
    this.postFxTextureB = this.core.device.createTexture({ size: this.textureSize, format: this.core.presentationFormat, usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC });
    this.orbitsTexture = this.core.device.createTexture({ size: this.textureSize, format: this.core.presentationFormat, usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });
    
    PostFXPass.clearTexture(this.core.device, this.postFxTextureA);
    PostFXPass.clearTexture(this.core.device, this.postFxTextureB);
    this.galaxyPass.onResize(this.textureSize, this.core);
  }

  private updateLoop = async (time: number) => {
    if (this.lastFrameTime === 0) this.lastFrameTime = time;
    const deltaTime = (time - this.lastFrameTime) / 1000.0;
    this.lastFrameTime = time;

    const clampedDeltaTime = Math.min(deltaTime, 1 / 20.0);
    this.lastDeltaTime = deltaTime;
    await this.authority.tick(clampedDeltaTime, { deltaX: this.input.deltaX, deltaY: this.input.deltaY, keys: this.input.keys });
    
    const systemState = await this.authority.query();
    this.lastSystemState = systemState;
    
    if (systemState.bodies.length !== this.scene.lastKnownBodyCount) {
        this.scene.recreateSpheresBuffer(systemState.bodies.length);
        this.computePass.recreatePipeline(systemState.bodies.length);
    }

    let bodiesToRender: Body[] = systemState.bodies;
    let renderScale: number;
    const camera = this.getCamera();

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
      this.cameraManager.update(this.state.cameraMode, { bodies: bodiesForMap, scale: renderScale, viewport: this.textureSize, vfov: 25.0, referenceFrame: this.state.referenceFrame });
      bodiesToRender = bodiesForMap;
      if (this.state.showOrbits) {
        // Recompute orbit geometry at correct scale
        this.orbitsPass.update(bodiesForMap, this.scene, this.core, renderScale);
        this.glyphsPass.update(
          this.orbitsPass.periapsisPoints,
          this.orbitsPass.apoapsisPoints,
          this.orbitsPass.ascendingNodePoints,
          this.orbitsPass.descendingNodePoints
        );
        this.soiPass.update(bodiesForMap, this.scene, renderScale);
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
    this.writeCameraBuffer(camera);

    this.scene.update(
      systemState,
      camera,
      bodiesToRender,
      renderScale,
      this.state.cameraMode === CameraMode.SYSTEM_MAP
    );
    // Do not update orbits in ship view; rings are hidden there

    // Handle click selection in System Map
    if (this.state.cameraMode === CameraMode.SYSTEM_MAP && this.input.clicked && this.textureSize) {
      const viewport = this.textureSize;
      const viewProj = camera.viewProjectionMatrix as unknown as number[];
      let closestBodyId: string | null = null;
      let closestDist2 = Number.POSITIVE_INFINITY;
      // Recompute the transformed bodies for hit testing, same as render path
      const fb = systemState.bodies.find(b => b.id === this.getCamera().focusBodyId);
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
        this.getCamera().focusBodyId = closestBodyId;
        this.getCamera().pendingFrame = true;
        if (this.ui) {
          const bodyName = systemState.bodies.find(b => b.id === closestBodyId)?.name;
          if (bodyName) {
            this.ui.setFocus(bodyName);
          }
        }
      }
    }

    this.render(renderScale);
    
    if (this.hud && this.textureSize) {
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
        this.hud.draw(hudBodies, camera, this.textureSize, this.state.cameraMode, this.state.playerShipId);
    }

    this.input.tick();
    requestAnimationFrame(this.updateLoop);
  }

  private writeCameraBuffer(camera: Camera) {
    const buffer = new Float32Array(64);
    buffer.set(camera.viewMatrix as unknown as number[], 0);
    buffer.set(camera.projectionMatrix as unknown as number[], 16);
    buffer.set(camera.viewProjectionMatrix as unknown as number[], 32);
    buffer[48] = camera.eye[0]; buffer[49] = camera.eye[1]; buffer[50] = camera.eye[2]; buffer[51] = 0.0;
    const dist = vec3.distance(camera.eye as unknown as number[], camera.look_at as unknown as number[]);
    buffer[52] = camera.forward[0]; buffer[53] = camera.forward[1]; buffer[54] = camera.forward[2]; buffer[55] = dist;
    buffer[56] = camera.right[0]; buffer[57] = camera.right[1]; buffer[58] = camera.right[2]; buffer[59] = 0.0;
    buffer[60] = camera.up[0]; buffer[61] = camera.up[1]; buffer[62] = camera.up[2]; buffer[63] = 0.0;
    this.core.device.queue.writeBuffer(this.scene.sharedCameraUniformBuffer, 0, buffer);
  }

  private render(systemScale: number) {
    if (!this.core.device || !this.textureSize) return;

    const theme = this.setTheme(this.currentThemeName, this.currentResponseName);
    if (!theme) return;

    const sourceTex = this.frameCount % 2 === 0 ? this.postFxTextureB : this.postFxTextureA;
    const destTex = this.frameCount % 2 === 0 ? this.postFxTextureA : this.postFxTextureB;

    const context: RenderContext = {
      core: this.core, scene: this.scene, camera: this.getCamera(), systemScale, textureSize: this.textureSize,
      sourceTexture: sourceTex, destinationTexture: destTex, mainSceneTexture: this.mainSceneTexture,
      orbitsTexture: this.orbitsTexture, themeUniformBuffer: this.themeUniformBuffer, lastDeltaTime: this.lastDeltaTime,
    };

    const encoder = this.core.device.createCommandEncoder();

    if (this.state.cameraMode === CameraMode.GALACTIC_MAP) {
        this.galaxyPass.run(encoder, context, this.galaxy.stars.length);
    } else if (this.state.cameraMode === CameraMode.SYSTEM_MAP) {
        if (this.state.showOrbits) {
            this.orbitsPass.run(encoder, context, theme);
        }
        this.mapPass.run(encoder, context, theme, this.state.showOrbits);
        if (this.state.showOrbits) {
          this.soiPass.run(encoder, context);
          this.glyphsPass.run(encoder, context);
        }
    } else { // This is now the Ops View (SHIP_RELATIVE)
        this.computePass.run(encoder, context);
        // Clear any stale orbits overlay when in ship view
        PostFXPass.clearTexture(this.core.device, this.orbitsTexture);
        this.postfxPass.run(encoder, context, theme, this.lastSystemState?.bodies ?? []);
    }

    this.core.device.queue.submit([encoder.finish()]);
    this.frameCount++;
  }
  
  public setTheme(themeName: string, responseName: string): Theme | undefined {
    this.currentThemeName = themeName;
    this.currentResponseName = responseName;
    const theme = themes[this.currentThemeName as keyof typeof themes];
    const response = spectralResponses[this.currentResponseName];
    if (!theme || !response) return undefined;

    const themeData = new Float32Array(20);
    themeData.set(theme.bg, 0); themeData.set(theme.fg, 4); themeData.set(theme.accent, 8);
    themeData.set(response, 12); themeData[16] = this.lastDeltaTime; themeData[17] = this.state.crtIntensity;
    this.core.device.queue.writeBuffer(this.themeUniformBuffer, 0, themeData);
    return theme;
  }

  public clearOrbitHistory(): void {
    this.orbitsPass.clearAll();
  }
}

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
import { mat4, vec4 } from 'gl-matrix';
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
import type { Body, Orbit, Ship, Theme, Vec3 } from '../shared/types';

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
  private frameCount = 0;
  
  private lastSystemState?: import('../shared/types').SystemState;

  // Orbit caching and draw state
  private orbitalElementsCache: Map<string, Orbit | null> = new Map();
  private lastFocusForCache: string | null = null;
  private currentOrbitsToDraw: (Orbit | null)[] = [];
  private currentTargetTimeSeconds: number = 0;
  // Note: prepared bodies/scale are derived each frame; keep locals instead of members

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

  private updateLoop = async () => {
    const deltaTime = 1 / 60.0;
    this.lastDeltaTime = deltaTime;
    await this.authority.tick(deltaTime, { deltaX: this.input.deltaX, deltaY: this.input.deltaY, keys: this.input.keys });
    
    const systemState = await this.authority.query();
    this.lastSystemState = systemState;
    
    if (systemState.bodies.length !== this.scene.lastKnownBodyCount) {
      this.scene.recreateSpheresBuffer(systemState.bodies.length);
      this.computePass.recreatePipeline(systemState.bodies.length);
      this.orbitsPass.resizeOrbitBuffers(this.core.device, systemState.bodies.length);
      // Reset current draw list to match size
      this.currentOrbitsToDraw = Array.from({ length: systemState.bodies.length }, () => null);
    }

    // Refresh scene hierarchy first (no draw data yet)
    this.scene.update(systemState, this.getCamera(), [], 1.0, false);

    let bodiesToRender: Body[] = systemState.bodies;
    let renderScale: number;

    if (this.state.cameraMode === CameraMode.SYSTEM_MAP) {
      // Invalidate cache when focus body changes
      const currentFocusId = this.getCamera().focusBodyId;
      if (currentFocusId !== this.lastFocusForCache) {
        this.orbitalElementsCache.clear();
        this.lastFocusForCache = currentFocusId;
      }

      const timestampSeconds = this.scene.lastSystemTimestamp / 1000.0;
      const hierarchy = this.scene.hierarchy;
      const bodies = systemState.bodies;
      const bodyMap = new Map(bodies.map(b => [b.id, b]));

      // Build orbits and predicted positions using analytic model
      const orbitsToDraw: (Orbit | null)[] = Array.from({ length: bodies.length }, () => null);
      const predictedBodies: Body[] = bodies.map((body, i) => {
        const parentId = hierarchy.get(body.id);
        const parent = parentId ? bodyMap.get(parentId) || null : null;
        if (!parent) {
          orbitsToDraw[i] = null;
          return body;
        }

        let orbit = this.orbitalElementsCache.get(body.id) ?? null;
        if (!orbit) {
          orbit = OrbitsPass.calculateOrbit(body, parent, timestampSeconds);
          this.orbitalElementsCache.set(body.id, orbit);
        }
        orbitsToDraw[i] = orbit;

        if (orbit) {
          const futureTime = timestampSeconds + this.state.timeOffset;
          const posInParentFrame = OrbitsPass.getPositionOnOrbit(orbit, futureTime);
          const newPosition: Vec3 = [
            parent.position[0] + posInParentFrame[0],
            parent.position[1] + posInParentFrame[1],
            parent.position[2] + posInParentFrame[2],
          ];
          return { ...body, position: newPosition };
        }
        return body;
      });

      // Apply reference frame to the predicted state
      let finalRenderBodies: Body[] = predictedBodies;
      const focusBody = finalRenderBodies.find(b => b.id === this.getCamera().focusBodyId);
      if (this.state.referenceFrame === ReferenceFrame.FOCUSED_BODY && focusBody) {
        finalRenderBodies = finalRenderBodies.map(body => ({
          ...body,
          position: [
            body.position[0] - focusBody.position[0],
            body.position[1] - focusBody.position[1],
            body.position[2] - focusBody.position[2],
          ] as Vec3,
        }));
      }

      // Render using final state
      renderScale = Scene.calculateRenderScale(systemState.bodies, this.getCamera().focusBodyId);
      this.cameraManager.update(this.state.cameraMode, { bodies: finalRenderBodies, scale: renderScale, viewport: this.textureSize, vfov: 25.0, referenceFrame: this.state.referenceFrame });
      bodiesToRender = finalRenderBodies;

      if (this.state.showOrbits) {
        // Refresh orbit geometry for current scale
        for (let i = 0; i < bodies.length; i++) {
          const orbit = orbitsToDraw[i];
          if (orbit) {
            this.orbitsPass.updateOrbitGeometry(this.core, i, orbit);
          }
        }
        this.currentOrbitsToDraw = orbitsToDraw;
        this.currentTargetTimeSeconds = timestampSeconds + this.state.timeOffset;

        this.glyphsPass.update(
          this.orbitsPass.periapsisPoints,
          this.orbitsPass.apoapsisPoints,
          this.orbitsPass.ascendingNodePoints,
          this.orbitsPass.descendingNodePoints
        );
        this.soiPass.update(finalRenderBodies, this.scene, renderScale);
      }
    } else if (this.state.cameraMode === CameraMode.SHIP_RELATIVE) {
      renderScale = 1.0;
      const playerShip = systemState.bodies.find(b => b.id === this.state.playerShipId) as Ship | undefined;
      this.cameraManager.update(this.state.cameraMode, { playerShip });
      bodiesToRender = systemState.bodies.filter(b => b.id !== this.state.playerShipId);
    } else {
      renderScale = 1.0;
      this.cameraManager.update(this.state.cameraMode, {});
    }

    this.scene.update(systemState, this.getCamera(), bodiesToRender, renderScale, this.state.cameraMode === CameraMode.SYSTEM_MAP);
    // Do not update orbits in ship view; rings are hidden there

    // Handle click selection in System Map (use final render bodies for hit testing)
    if (this.state.cameraMode === CameraMode.SYSTEM_MAP && this.input.clicked && this.textureSize) {
      const viewport = this.textureSize;
      const viewProj = mat4.multiply(mat4.create(), this.getCamera().projectionMatrix as unknown as number[], this.getCamera().viewMatrix as unknown as number[]);
      let closestBodyId: string | null = null;
      let closestDist2 = Number.POSITIVE_INFINITY;
      for (const b of bodiesToRender) {
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
          const bodyName = bodiesToRender.find(b => b.id === closestBodyId)?.name;
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
          const scaledFinalBodies = bodiesToRender.map(b => ({
            ...b,
            position: [ b.position[0] * renderScale, b.position[1] * renderScale, b.position[2] * renderScale ] as Vec3,
          }));
          hudBodies = scaledFinalBodies;
        }
        this.hud.draw(hudBodies, this.getCamera(), this.textureSize, this.state.cameraMode, this.state.playerShipId);
    }

    this.input.tick();
    requestAnimationFrame(this.updateLoop);
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
          this.orbitsPass.run(encoder, context, theme, this.currentOrbitsToDraw, this.currentTargetTimeSeconds);
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

  public clearOrbitCache(): void {
    this.orbitalElementsCache.clear();
  }
}

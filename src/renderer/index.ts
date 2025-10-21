import type { Authority } from '../authority/authority';
import type { InputManager } from '../input';
import { AppState, CameraMode } from '../state';
import { CameraManager } from '../camera/manager';
import type { Camera } from '../camera';
import { HUDManager } from '../hud';
import { Galaxy } from '../galaxy';
import { WebGPUCore } from './core';
import { Scene } from './scene';
import type { RenderContext } from './types';
import { ComputePass } from './passes/computePass';
import { GalaxyPass } from './passes/galaxyPass';
import { OrbitsPass } from './passes/orbitsPass';
import { PostFXPass } from './passes/postfxPass';
import { themes } from '../theme';
import { spectralResponses } from '../spectral';
import type { Body, Ship, Theme } from '../shared/types';

export class Renderer {
  public readonly authority: Authority;
  private readonly canvas: HTMLCanvasElement;
  private readonly state: AppState;
  private readonly input: InputManager;
  
  private core!: WebGPUCore;
  private scene!: Scene;
  private cameraManager: CameraManager;
  private hud!: HUDManager;
  private galaxy: Galaxy;
  
  private computePass!: ComputePass;
  private galaxyPass!: GalaxyPass;
  private orbitsPass!: OrbitsPass;
  private postfxPass!: PostFXPass;

  private textureSize!: { width: number, height: number };
  private mainSceneTexture!: GPUTexture;
  private postFxTextureA!: GPUTexture;
  private postFxTextureB!: GPUTexture;
  private orbitsTexture!: GPUTexture;
  
  private themeUniformBuffer!: GPUBuffer;
  private currentThemeName: string = 'amber';
  private currentResponseName: string = 'Visible (Y)';
  private lastDeltaTime: number = 1 / 60;
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
    
    this.computePass.initialize(this.core, this.scene);
    this.galaxyPass.initialize(this.core, this.scene);
    this.orbitsPass.initialize(this.core, this.scene);
    this.postfxPass.initialize(this.core, this.scene);

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
    this.input.tick();
    
    const systemState = await this.authority.query();
    this.lastSystemState = systemState;
    
    if (systemState.bodies.length !== this.scene.lastKnownBodyCount) {
        this.scene.recreateSpheresBuffer(systemState.bodies.length);
        this.computePass.recreatePipeline(systemState.bodies.length);
    }

    let bodiesToRender: Body[] = systemState.bodies;
    let renderScale: number;

    if (this.state.cameraMode === CameraMode.SYSTEM_ORBITAL) {
      renderScale = Scene.calculateRenderScale(systemState.bodies);
      this.cameraManager.update(this.state.cameraMode, { bodies: systemState.bodies, scale: renderScale, viewport: this.textureSize, vfov: 25.0 });
    } else if (this.state.cameraMode === CameraMode.SHIP_RELATIVE) {
      renderScale = 1.0;
      const playerShip = systemState.bodies.find(b => b.id === this.state.playerShipId) as Ship | undefined;
      this.cameraManager.update(this.state.cameraMode, { playerShip });
      bodiesToRender = systemState.bodies.filter(b => b.id !== this.state.playerShipId);
    } else {
      renderScale = 1.0;
      this.cameraManager.update(this.state.cameraMode, {});
    }

    this.scene.update(systemState, this.getCamera(), bodiesToRender, renderScale);
    if (this.state.showOrbits) {
      this.orbitsPass.update(systemState.bodies, this.scene, this.core, renderScale);
    }

    this.render(renderScale);
    
    if (this.hud && this.textureSize) {
        const scaledBodies = systemState.bodies.map(b => ({
            ...b, position: [ b.position[0] * renderScale, b.position[1] * renderScale, b.position[2] * renderScale ] as [number,number,number]
        }));
        this.hud.draw(this.state.cameraMode === CameraMode.SYSTEM_ORBITAL ? scaledBodies : systemState.bodies, this.getCamera(), this.textureSize, this.state.cameraMode, this.state.playerShipId);
    }

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
    } else {
      this.computePass.run(encoder, context);
      if (this.state.showOrbits) {
        this.orbitsPass.run(encoder, context, theme);
      }
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
    themeData.set(response, 12); themeData[16] = this.lastDeltaTime;
    this.core.device.queue.writeBuffer(this.themeUniformBuffer, 0, themeData);
    return theme;
  }

  public clearOrbitHistory(): void {
    this.orbitsPass.clearAll();
  }
}

import type { Authority } from '../authority';
import { AppState, CameraMode } from '../state';
import type { Camera } from '../camera';
import { WebGPUCore } from './core';
import { Scene } from './scene';
import type { RenderContext } from './types';
import { vec3 } from 'gl-matrix';
import type { UI } from '../ui';
import { PostFXPass } from './passes/postfxPass';
import { themes } from '../theme';
import { spectralResponses } from '../spectral';
import type { Theme, FrameData } from '@shared/types';
import { HUDManager } from '../hud';
import type { SystemState } from '@shared/types';
import type { IRenderPipeline } from './pipelines/base';

export class Renderer {
  private readonly canvas: HTMLCanvasElement;
  private readonly state: AppState;
  public ui: UI | null = null;
  private readonly hud: HUDManager;
  
  private core!: WebGPUCore;
  private scene!: Scene;
  

  private textureSize!: { width: number, height: number };
  private mainSceneTexture!: GPUTexture;
  private postFxTextureA!: GPUTexture;
  private postFxTextureB!: GPUTexture;
  private orbitsTexture!: GPUTexture;

  private themeUniformBuffer!: GPUBuffer;
  private sceneUniformBuffer!: GPUBuffer;
  private pipelines!: Record<CameraMode, IRenderPipeline>;
  private currentThemeName: string = 'white';
  private currentResponseName: string = 'Full Color';
  private lastDeltaTime: number = 1 / 60;
  private lastFrameTime: number = 0; // TODO: consider removing if unused
  private frameCount = 0;
  
  private lastSystemState?: SystemState; // TODO: consider removing if unused

  constructor(
    canvas: HTMLCanvasElement,
    state: AppState,
    hud: HUDManager,
    pipelines: Record<CameraMode, IRenderPipeline>
  ) {
    this.canvas = canvas;
    this.state = state;
    this.hud = hud;
    this.pipelines = pipelines;
  }

  public async initialize(authority: Authority): Promise<void> {
    this.core = await WebGPUCore.create(this.canvas);
    this.scene = new Scene(this.core);

    const initialSystemState = await authority.query();
    this.lastSystemState = initialSystemState;
    this.scene.initialize(initialSystemState, []);
    
    this.themeUniformBuffer = this.core.device.createBuffer({
      size: 80,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.sceneUniformBuffer = this.core.device.createBuffer({
      size: 48, // three vec4<f32> (added tier_scale_and_pad)
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    for (const pipeline of Object.values(this.pipelines)) {
      await pipeline.initialize(this.core, this.scene);
    }

    this.handleResize();
    
    new ResizeObserver(this.handleResize).observe(this.canvas);
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

    this.mainSceneTexture = this.core.device.createTexture({ size: this.textureSize, format: 'rgba16float', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT });
    this.postFxTextureA = this.core.device.createTexture({ size: this.textureSize, format: 'rgba16float', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC });
    this.postFxTextureB = this.core.device.createTexture({ size: this.textureSize, format: 'rgba16float', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC });
    this.orbitsTexture = this.core.device.createTexture({ size: this.textureSize, format: this.core.presentationFormat, usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });

    PostFXPass.clearTexture(this.core.device, this.postFxTextureA);
    PostFXPass.clearTexture(this.core.device, this.postFxTextureB);
    
    for (const pipeline of Object.values(this.pipelines)) {
      pipeline.onResize(this.textureSize, this.core);
    }
  }

  

  public writeCameraBuffer(camera: Camera) {
    const buffer = new Float32Array(80);
    buffer.set(camera.viewMatrix as unknown as number[], 0);
    buffer.set(camera.projectionMatrix as unknown as number[], 16);
    buffer.set(camera.viewProjectionMatrix as unknown as number[], 32);
    buffer[48] = camera.eye[0]; buffer[49] = camera.eye[1]; buffer[50] = camera.eye[2]; buffer[51] = 0.0;
    const dist = vec3.distance(camera.eye as unknown as number[], camera.look_at as unknown as number[]);
    buffer[52] = camera.forward[0]; buffer[53] = camera.forward[1]; buffer[54] = camera.forward[2]; buffer[55] = dist;
    buffer[56] = camera.right[0]; buffer[57] = camera.right[1]; buffer[58] = camera.right[2]; buffer[59] = 0.0;
    buffer[60] = camera.up[0]; buffer[61] = camera.up[1]; buffer[62] = camera.up[2]; buffer[63] = 0.0;
    // projection_constants (vec4): .x = lod_constant
    const vfov_rad = camera.vfov * (Math.PI / 180.0);
    const lod_constant = this.textureSize.height / vfov_rad;
    buffer[64] = lod_constant;
    this.core.device.queue.writeBuffer(this.scene.sharedCameraUniformBuffer, 0, buffer);
  }

  public getCore(): WebGPUCore { return this.core; }
  public getScene(): Scene { return this.scene; }
  public getTextureSize(): { width: number, height: number } { return this.textureSize; }
  public getCanvas(): HTMLCanvasElement { return this.canvas; }

  public prepare(frameData: FrameData) {
    const pipeline = this.pipelines[frameData.cameraMode as unknown as keyof typeof this.pipelines];
    pipeline?.prepare(frameData, this);
  }

  public render(frameData: FrameData) {
    if (!this.core.device || !this.textureSize) return;
    const { camera, systemScale, deltaTime } = frameData;
    this.lastDeltaTime = deltaTime;
    const theme = this.setTheme(this.currentThemeName, this.currentResponseName);
    if (!theme) return;

    const sourceTex = this.frameCount % 2 === 0 ? this.postFxTextureB : this.postFxTextureA;
    const destTex = this.frameCount % 2 === 0 ? this.postFxTextureA : this.postFxTextureB;

    const context: RenderContext = {
      core: this.core,
      scene: this.scene,
      camera: camera,
      systemScale,
      textureSize: this.textureSize,
      sourceTexture: sourceTex,
      destinationTexture: destTex,
      mainSceneTexture: this.mainSceneTexture,
      orbitsTexture: this.orbitsTexture,
      themeUniformBuffer: this.themeUniformBuffer,
      sceneUniformBuffer: this.sceneUniformBuffer,
      lastDeltaTime: deltaTime,
    };

    const encoder = this.core.device.createCommandEncoder();

    const pipeline = this.pipelines[frameData.cameraMode as unknown as keyof typeof this.pipelines];
    pipeline?.render(encoder, context, frameData, theme);

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
    themeData.set(response, 12); 
    themeData[16] = this.lastDeltaTime; 
    themeData[17] = this.state.crtIntensity;
    this.core.device.queue.writeBuffer(this.themeUniformBuffer, 0, themeData);
    return theme;
  }

  public clearOrbitHistory(): void {
    const systemMap = this.pipelines[CameraMode.SYSTEM_MAP as unknown as keyof typeof this.pipelines] as any;
    systemMap?.clearOrbitHistory?.();
  }
}

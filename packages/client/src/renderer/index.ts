import type { Authority } from '../authority';
import { AppState, CameraMode } from '../state';
import type { Camera } from '../camera/camera';
import { WebGPUCore } from './core';
import { Scene } from './scene';
import type { RenderContext } from './types';
import { vec3 } from 'gl-matrix';
import type { UI } from '../ui';
import { PostFXPass } from './postfx/postfxPass';
import { HydraulicsPass } from './hydraulics/hydraulicsPass';
import { TerrainBakePass } from './hydraulics/terrainBakePass';
import { themes } from '../theme';
import { spectralResponses } from '../spectral';
import type { Theme } from '@shared/types';
import { HUDManager } from '../hud';
import type { SystemState } from '@shared/types';
import type { IRenderPipeline } from './pipelines/base';
import { WasmGridBridge } from '../authority/wasmGrid';
 
import type { RenderPayload } from '@client/orchestration/types';
import { SceneDataProcessor } from './data/sceneDataProcessor';
import type { Vec3 } from '@shared/types';

export class Renderer {
  private readonly canvas: HTMLCanvasElement;
  private readonly state: AppState;
  public ui: UI | null = null;
  private readonly hud: HUDManager;
  
  private core!: WebGPUCore;
  private scene!: Scene;
  private hydraulicsPass!: HydraulicsPass;
  private terrainBakePass!: TerrainBakePass;
  

  private textureSize!: { width: number, height: number };
  private mainSceneTexture!: GPUTexture;
  private postFxTextureA!: GPUTexture;
  private postFxTextureB!: GPUTexture;
  private orbitsTexture!: GPUTexture;
  private terrainHeightTexture!: GPUTexture;
  private waterStateTextureA!: GPUTexture;
  private waterStateTextureB!: GPUTexture;

  private themeUniformBuffer!: GPUBuffer;
  private sceneUniformBuffer!: GPUBuffer;
  private pipelines!: Record<CameraMode, IRenderPipeline>;
  private currentThemeName: string = 'white';
  private currentResponseName: string = 'Full Color';
  private lastDeltaTime: number = 1 / 60;
  private lastFrameTime: number = 0; // TODO: consider removing if unused
  private frameCount = 0;
  
  private lastSystemState?: SystemState; // TODO: consider removing if unused
  private sceneDataProcessor: SceneDataProcessor;

  // Coarse grid GPU buffers
  private gridVertexBuffer!: GPUBuffer;
  private gridElevationBuffer!: GPUBuffer;
  private gridIndexBuffer!: GPUBuffer;

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
    this.sceneDataProcessor = new SceneDataProcessor();
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
      size: 48, // three vec4<f32>
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    this.hydraulicsPass = new HydraulicsPass();
    await this.hydraulicsPass.initialize(this.core, this.scene);
    this.terrainBakePass = new TerrainBakePass();
    await this.terrainBakePass.initialize(this.core, this.scene);

    // Initialize coarse grid and upload to GPU once
    await this.initializeGrid(3);

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
      const dpr = (globalThis as any).devicePixelRatio || 1;
      this.hud.resize(width, height, dpr);
    }
    this.textureSize = { width, height };
    this.core.context.configure({ device: this.core.device, format: this.core.presentationFormat });

    [this.mainSceneTexture, this.postFxTextureA, this.postFxTextureB, this.orbitsTexture, this.terrainHeightTexture, this.waterStateTextureA, this.waterStateTextureB]
        .forEach(tex => tex?.destroy());

    this.mainSceneTexture = this.core.device.createTexture({ size: this.textureSize, format: 'rgba16float', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT });
    this.postFxTextureA = this.core.device.createTexture({ size: this.textureSize, format: 'rgba16float', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC });
    this.postFxTextureB = this.core.device.createTexture({ size: this.textureSize, format: 'rgba16float', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC });
    this.orbitsTexture = this.core.device.createTexture({ size: this.textureSize, format: this.core.presentationFormat, usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });

    // Simulation grid textures are fixed size, independent of screen resolution
    const simGridSize = { width: 1024, height: 1024, depthOrArrayLayers: 6 } as const;
    const simTextureUsage = GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC;
    this.terrainHeightTexture = this.core.device.createTexture({ size: simGridSize, format: 'r32float', usage: simTextureUsage });
    this.waterStateTextureA = this.core.device.createTexture({ size: simGridSize, format: 'rgba16float', usage: simTextureUsage });
    this.waterStateTextureB = this.core.device.createTexture({ size: simGridSize, format: 'rgba16float', usage: simTextureUsage });

    // Bake terrain and initialize water based on sea level
    const tempContext: RenderContext = {
      core: this.core,
      scene: this.scene,
      camera: undefined as unknown as Camera, // not used by bake/init
      systemScale: 1.0,
      textureSize: this.textureSize,
      sourceTexture: this.postFxTextureA,
      destinationTexture: this.postFxTextureB,
      mainSceneTexture: this.mainSceneTexture,
      orbitsTexture: this.orbitsTexture,
      themeUniformBuffer: this.themeUniformBuffer,
      sceneUniformBuffer: this.sceneUniformBuffer,
      lastDeltaTime: 0,
      gridVertexBuffer: this.gridVertexBuffer,
      gridElevationBuffer: this.gridElevationBuffer,
      gridIndexBuffer: this.gridIndexBuffer,
      terrainHeight: this.terrainHeightTexture,
      waterRead: this.waterStateTextureA,
      waterWrite: this.waterStateTextureB,
    } as unknown as RenderContext;
    this.terrainBakePass.run(tempContext, this.terrainHeightTexture);
    this.hydraulicsPass.initializeState(this.waterStateTextureA, tempContext, this.terrainHeightTexture);
    this.hydraulicsPass.initializeState(this.waterStateTextureB, tempContext, this.terrainHeightTexture);

    PostFXPass.clearTexture(this.core.device, this.postFxTextureA);
    PostFXPass.clearTexture(this.core.device, this.postFxTextureB);
    
    for (const pipeline of Object.values(this.pipelines)) {
      pipeline.onResize(this.textureSize, this.core);
    }
  }

  private async initializeGrid(size: number): Promise<void> {
    const bridge = await WasmGridBridge.create();
    bridge.createGrid(size);
    const vertices = bridge.getGridVertexBuffer(); // Float32Array of xyz triplets
    const elevations = bridge.getGridElevationBuffer(); // Float32Array per-vertex
    const indices = bridge.getGridIndexBuffer(); // Uint32Array triangle indices

    const device = this.core.device;
    this.gridVertexBuffer = device.createBuffer({
      label: 'Grid Vertex Buffer',
      size: Math.max(4, vertices.byteLength),
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.gridElevationBuffer = device.createBuffer({
      label: 'Grid Elevation Buffer',
      size: Math.max(4, elevations.byteLength),
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.gridIndexBuffer = device.createBuffer({
      label: 'Grid Index Buffer',
      size: Math.max(4, indices.byteLength),
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

    device.queue.writeBuffer(this.gridVertexBuffer, 0, vertices.buffer, vertices.byteOffset, vertices.byteLength);
    device.queue.writeBuffer(this.gridElevationBuffer, 0, elevations.buffer, elevations.byteOffset, elevations.byteLength);
    device.queue.writeBuffer(this.gridIndexBuffer, 0, indices.buffer, indices.byteOffset, indices.byteLength);
  }

  

  public writeCameraBuffer(camera: Camera, worldCameraEye: Vec3) {
    // The new buffer has 88 floats (352 bytes)
    const buffer = new Float32Array(88);
    buffer.set(camera.viewMatrix as unknown as number[], 0);
    buffer.set(camera.projectionMatrix as unknown as number[], 16);
    buffer.set(camera.viewProjectionMatrix as unknown as number[], 32);

    // Split the f64 world position and write to the new slots
    const [eyeX_h, eyeX_l] = this.sceneDataProcessor.splitDouble(worldCameraEye[0] as number);
    const [eyeY_h, eyeY_l] = this.sceneDataProcessor.splitDouble(worldCameraEye[1] as number);
    const [eyeZ_h, eyeZ_l] = this.sceneDataProcessor.splitDouble(worldCameraEye[2] as number);
    buffer[48] = eyeX_h; buffer[49] = eyeY_h; buffer[50] = eyeZ_h; // eye_pos_high
    buffer[52] = eyeX_l; buffer[53] = eyeY_l; buffer[54] = eyeZ_l; // eye_pos_low

    // The old `eye` uniform is now the camera-relative origin, always (0,0,0)
    buffer[56] = 0.0; buffer[57] = 0.0; buffer[58] = 0.0;

    const dist = vec3.distance(camera.eye as unknown as number[], camera.look_at as unknown as number[]);
    buffer[60] = camera.forward[0]; buffer[61] = camera.forward[1]; buffer[62] = camera.forward[2]; buffer[63] = dist;
    buffer[64] = camera.right[0]; buffer[65] = camera.right[1]; buffer[66] = camera.right[2];
    buffer[68] = camera.up[0]; buffer[69] = camera.up[1]; buffer[70] = camera.up[2];
    
    const vfov_rad = camera.vfov * (Math.PI / 180.0);
    const lod_constant = this.textureSize.height / vfov_rad;
    buffer[72] = lod_constant;

    this.core.device.queue.writeBuffer(this.scene.sharedCameraUniformBuffer, 0, buffer);
  }

  public getCore(): WebGPUCore { return this.core; }
  public getScene(): Scene { return this.scene; }
  public getTextureSize(): { width: number, height: number } { return this.textureSize; }
  public getCanvas(): HTMLCanvasElement { return this.canvas; }

  public prepare(frameData: RenderPayload) {
    const pipeline = this.pipelines[frameData.cameraMode as unknown as keyof typeof this.pipelines];
    pipeline?.prepare(frameData, this);
  }

  public render(frameData: RenderPayload) {
    if (!this.core.device || !this.textureSize) return;
    const { camera, deltaTime } = frameData;
    const systemScaleValue = 1.0; // System Map removed
    this.lastDeltaTime = deltaTime;
    const theme = this.setTheme(this.currentThemeName, this.currentResponseName);
    if (!theme) return;

    // Ensure theme buffer exists (defensive)
    if (!this.themeUniformBuffer) {
      this.themeUniformBuffer = this.core.device.createBuffer({ size: 80, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    }

    const sourceTex = this.frameCount % 2 === 0 ? this.postFxTextureB : this.postFxTextureA;
    const destTex = this.frameCount % 2 === 0 ? this.postFxTextureA : this.postFxTextureB;

    // Ping-pong simulation textures
    const waterReadTex = this.frameCount % 2 === 0 ? this.waterStateTextureA : this.waterStateTextureB;
    const waterWriteTex = this.frameCount % 2 === 0 ? this.waterStateTextureB : this.waterStateTextureA;

    const context: RenderContext = {
      core: this.core,
      scene: this.scene,
      camera: camera,
      systemScale: systemScaleValue,
      textureSize: this.textureSize,
      sourceTexture: sourceTex,
      destinationTexture: destTex,
      mainSceneTexture: this.mainSceneTexture,
      orbitsTexture: this.orbitsTexture,
      themeUniformBuffer: this.themeUniformBuffer,
      sceneUniformBuffer: this.sceneUniformBuffer,
      lastDeltaTime: deltaTime,
      gridVertexBuffer: this.gridVertexBuffer,
      gridElevationBuffer: this.gridElevationBuffer,
      gridIndexBuffer: this.gridIndexBuffer,
      waterRead: waterReadTex,
      waterWrite: waterWriteTex,
    };

    const encoder = this.core.device.createCommandEncoder();

    // Pass 1: Run hydraulics simulation
    this.hydraulicsPass.run(encoder, context, this.terrainHeightTexture);

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
}

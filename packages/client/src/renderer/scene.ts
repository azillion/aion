import type { Body, SystemState, Star } from '@shared/types';
import type { WebGPUCore } from './core';
import type { Camera } from '../camera/camera';

const FLOATS_PER_SPHERE = 32; // This was correctly updated in Phase 1
const FLOATS_PER_STAR = 8;
const VISUAL_SETTINGS = { systemViewSize: 20.0 } as const;
const CAMERA_UNIFORM_BUFFER_SIZE = 352; // 88 floats * 4 bytes/float. Matches Float32Array(88) in renderer/index.ts

// Determine parent-child relationships by strongest gravitational influence
function buildSystemHierarchy(bodies: Body[]): Map<string, string | null> {
  const parentMap = new Map<string, string | null>();
  for (let i = 0; i < bodies.length; i++) {
    const body = bodies[i];
    let bestParentId: string | null = null;
    let maxForce = -Infinity;
    for (let j = 0; j < bodies.length; j++) {
      if (i === j) continue;
      const potentialParent = bodies[j];
      const dx = potentialParent.position[0] - body.position[0];
      const dy = potentialParent.position[1] - body.position[1];
      const dz = potentialParent.position[2] - body.position[2];
      const distSq = dx * dx + dy * dy + dz * dz;
      if (distSq <= 0) continue;
      const force = potentialParent.mass / distSq;
      if (force > maxForce) {
        maxForce = force;
        bestParentId = potentialParent.id;
      }
    }
    parentMap.set(body.id, bestParentId);
  }
  return parentMap;
}

export class Scene {
  private device: GPUDevice;

  // Scene object buffers
  public sceneObjectsBuffer!: GPUBuffer;
  public mapSpheresBuffer!: GPUBuffer;
  public shadowCasterBuffer!: GPUBuffer;
  public shadowCasterCountBuffer!: GPUBuffer;
  public starBuffer!: GPUBuffer;

  // A SINGLE, shared buffer for all camera data.
  public sharedCameraUniformBuffer!: GPUBuffer;

  public lastKnownBodyCount: number = 0;
  public sceneObjectCount: number = 0;
  
  public get hierarchy(): Map<string, string | null> { return this._hierarchy; }
  private _hierarchy: Map<string, string | null> = new Map();

  constructor(core: WebGPUCore) {
    this.device = core.device;
    this.sharedCameraUniformBuffer = this.device.createBuffer({
      label: 'Shared Camera Uniform Buffer',
      size: CAMERA_UNIFORM_BUFFER_SIZE,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.shadowCasterBuffer = this.device.createBuffer({
      label: 'Shadow Caster Buffer',
      // Allocate room for up to 32 casters
      size: 32 * FLOATS_PER_SPHERE * 4,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.shadowCasterCountBuffer = this.device.createBuffer({
      label: 'Shadow Caster Count Buffer',
      size: 4,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }

  public initialize(systemState: SystemState, stars: Star[]) {
    this.lastKnownBodyCount = systemState.bodies.length;
    this.initializeSceneObjectBuffer(this.lastKnownBodyCount);

    const starData = this.serializeStars(stars);
    this.starBuffer = this.device.createBuffer({
      size: starData.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(this.starBuffer, 0, starData);
  }

  public updateSceneObjects(data: Float32Array, count: number) {
    if (count > this.lastKnownBodyCount) {
      this.initializeSceneObjectBuffer(count);
    }
    this.sceneObjectCount = count;
    this.device.queue.writeBuffer(this.sceneObjectsBuffer, 0, data.buffer, data.byteOffset, data.byteLength);
  }

  public updateShadowCasters(data: Float32Array, count: number) {
    const capacityFloats = this.shadowCasterBuffer.size / 4;
    const scratch = new Float32Array(capacityFloats);
    scratch.fill(0);
    scratch.set(data);
    this.device.queue.writeBuffer(this.shadowCasterBuffer, 0, scratch.buffer, 0, scratch.byteLength);
    this.device.queue.writeBuffer(this.shadowCasterCountBuffer, 0, new Uint32Array([count]));
  }

  public initializeSceneObjectBuffer(initialBodyCount: number) {
    const size = Math.ceil(initialBodyCount * 1.5) * FLOATS_PER_SPHERE * 4;

    if (this.sceneObjectsBuffer) this.sceneObjectsBuffer.destroy();
    this.sceneObjectsBuffer = this.device.createBuffer({ label: `Scene Objects Buffer`, size, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
    this.lastKnownBodyCount = initialBodyCount;
  }

  // Map view buffer management
  public recreateMapSpheresBuffer(numBodies: number) {
    if (this.mapSpheresBuffer) {
      this.mapSpheresBuffer.destroy();
    }
    const size = Math.max(1, numBodies * FLOATS_PER_SPHERE * 4);
    this.mapSpheresBuffer = this.device.createBuffer({
      label: `Map Spheres Buffer (${numBodies} bodies)`,
      size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.lastKnownBodyCount = numBodies;
  }

  public updateMapBuffer(systemState: SystemState, _camera: Camera, bodiesToRender: Body[], _renderScale: number, _useWorldSpace: boolean = false) {
    if (systemState.bodies.length !== this.lastKnownBodyCount || !this.mapSpheresBuffer) {
      this.recreateMapSpheresBuffer(systemState.bodies.length);
    }
    this._hierarchy = buildSystemHierarchy(systemState.bodies);
    const sphereData = this.serializeSystemState(bodiesToRender);
    this.device.queue.writeBuffer(this.mapSpheresBuffer, 0, sphereData);
  }

  private serializeSystemState(bodies: Body[]) {
    const sphereData = new Float32Array(bodies.length * FLOATS_PER_SPHERE);
    bodies.forEach((body, i) => {
      const f_base = i * FLOATS_PER_SPHERE;
      // In an f64 architecture, we would split the f64 position here.
      // However, the map view still uses a scaled f32 system.
      // This serialization logic is now specific to the map view.
      sphereData[f_base + 0] = body.position[0]; // pos_high.x
      sphereData[f_base + 1] = body.position[1]; // pos_high.y
      sphereData[f_base + 2] = body.position[2]; // pos_high.z
      sphereData[f_base + 3] = body.radius; // radius

      // pos_low is zero for the map view

      // Material
      sphereData.set(body.albedo, f_base + 8); // offset 8
      sphereData[f_base + 11] = body.terrain && (body.terrain as any).atmosphere ? 1.0 : 0.0;

      const emissive = body.emissive ?? [0, 0, 0];
      sphereData.set(emissive, f_base + 12); // offset 12
      sphereData[f_base + 15] = body.terrain ? 1.0 : 0.0;
      
      // Maintain struct alignment: ref_idx_opacity_pad vec4 slot
      // Even if unused, fill with sane defaults to match WGSL layout.
      sphereData[f_base + 16] = 1.0; // ref_idx (offset 16)
      sphereData[f_base + 17] = 1.0; // opacity

      // Planet data
      if (body.terrain) {
        sphereData[f_base + 20] = body.terrain.radius;      // base_radius (offset 20)
        sphereData[f_base + 21] = body.terrain.seaLevel;    // sea_level
        sphereData[f_base + 22] = body.terrain.maxHeight;   // max_height
        sphereData[f_base + 23] = body.terrain.noiseSeed;   // seed
      } else {
        // This also correctly zeroes the terrain_params block
        sphereData.fill(0.0, f_base + 20, f_base + 24);
      }
    });
    return sphereData;
  }
  
  public static calculateRenderScale(bodies: Body[], focusBodyId: string | null): number {
    const focusBody = bodies.find(b => b.id === focusBodyId);
    const hierarchy = buildSystemHierarchy(bodies);
    let systemRadius = 0;

    if (focusBody && focusBody.name !== 'Sun') {
      const children = bodies.filter(b => hierarchy.get(b.id) === focusBody.id);
      if (children.length > 0) {
        systemRadius = children.reduce((maxDist, b) => {
          const dist = Math.hypot(
            b.position[0] - focusBody.position[0],
            b.position[1] - focusBody.position[1],
            b.position[2] - focusBody.position[2]
          );
          return Math.max(maxDist, dist);
        }, 0);
      } else {
        systemRadius = focusBody.radius * 500;
      }
    } else {
      systemRadius = bodies.reduce((maxDist, b) => {
        const dist = Math.hypot(b.position[0], b.position[1], b.position[2]);
        return Math.max(maxDist, dist);
      }, 0);
    }
    
    return systemRadius > 0 ? VISUAL_SETTINGS.systemViewSize / systemRadius : 1.0;
  }

  private serializeStars(stars: Star[]) {
    const data = new Float32Array(stars.length * FLOATS_PER_STAR);
    for (let i = 0; i < stars.length; i++) {
      const base = i * FLOATS_PER_STAR;
      const s = stars[i];
      data.set(s.position, base + 0);
      data.set(s.color, base + 4);
      data[base + 7] = s.size;
    }
    return data;
  }
}


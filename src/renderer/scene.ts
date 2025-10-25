import type { Body, SystemState, Star } from '../shared/types';
import type { WebGPUCore } from './core';
import type { Camera } from '../camera';

const FLOATS_PER_SPHERE = 24;
const FLOATS_PER_STAR = 8;
const VISUAL_SETTINGS = { systemViewSize: 20.0 } as const;

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
  public nearTierBuffer!: GPUBuffer;
  public midTierBuffer!: GPUBuffer;
  public farTierBuffer!: GPUBuffer;
  public mapSpheresBuffer!: GPUBuffer;
  public starBuffer!: GPUBuffer;

  // A SINGLE, shared buffer for all camera data.
  public sharedCameraUniformBuffer!: GPUBuffer;

  public lastKnownBodyCount: number = 0;
  public nearCount: number = 0;
  public midCount: number = 0;
  public farCount: number = 0;
  // Reusable zero-padded CPU buffers matching GPU tier buffer sizes to avoid per-frame allocations
  private nearZeroPadded?: Float32Array;
  private midZeroPadded?: Float32Array;
  private farZeroPadded?: Float32Array;
  
  public get hierarchy(): Map<string, string | null> { return this._hierarchy; }
  private _hierarchy: Map<string, string | null> = new Map();

  constructor(core: WebGPUCore) {
    this.device = core.device;
    this.sharedCameraUniformBuffer = this.device.createBuffer({
      label: 'Shared Camera Uniform Buffer',
      size: 320,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }

  public initialize(systemState: SystemState, stars: Star[]) {
    this.lastKnownBodyCount = systemState.bodies.length;
    this.initializeTierBuffers(this.lastKnownBodyCount);

    const starData = this.serializeStars(stars);
    this.starBuffer = this.device.createBuffer({
      size: starData.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(this.starBuffer, 0, starData);
  }

  public updateTiers(nearBodies: Body[], midBodies: Body[], farBodies: Body[]) {
    const totalBodyCount = nearBodies.length + midBodies.length + farBodies.length;
    if (totalBodyCount > this.lastKnownBodyCount) {
      this.initializeTierBuffers(totalBodyCount);
    }
    this.nearCount = nearBodies.length;
    this.midCount = midBodies.length;
    this.farCount = farBodies.length;

    // Write full zero-padded buffers each frame using reusable scratch arrays
    {
      const data = this.serializeSystemState(nearBodies);
      if (!this.nearZeroPadded || this.nearZeroPadded.length !== (this.nearTierBuffer.size / 4)) {
        this.nearZeroPadded = new Float32Array(this.nearTierBuffer.size / 4);
      }
      this.nearZeroPadded.fill(0);
      this.nearZeroPadded.set(data);
      this.device.queue.writeBuffer(this.nearTierBuffer, 0, this.nearZeroPadded.buffer, 0, this.nearZeroPadded.byteLength);
    }
    {
      const data = this.serializeSystemState(midBodies);
      if (!this.midZeroPadded || this.midZeroPadded.length !== (this.midTierBuffer.size / 4)) {
        this.midZeroPadded = new Float32Array(this.midTierBuffer.size / 4);
      }
      this.midZeroPadded.fill(0);
      this.midZeroPadded.set(data);
      this.device.queue.writeBuffer(this.midTierBuffer, 0, this.midZeroPadded.buffer, 0, this.midZeroPadded.byteLength);
    }
    {
      const data = this.serializeSystemState(farBodies);
      if (!this.farZeroPadded || this.farZeroPadded.length !== (this.farTierBuffer.size / 4)) {
        this.farZeroPadded = new Float32Array(this.farTierBuffer.size / 4);
      }
      this.farZeroPadded.fill(0);
      this.farZeroPadded.set(data);
      this.device.queue.writeBuffer(this.farTierBuffer, 0, this.farZeroPadded.buffer, 0, this.farZeroPadded.byteLength);
    }
  }

  public initializeTierBuffers(initialBodyCount: number) {
    const size = Math.max(1, initialBodyCount * FLOATS_PER_SPHERE * 4);

    if (this.nearTierBuffer) this.nearTierBuffer.destroy();
    this.nearTierBuffer = this.device.createBuffer({
      label: `Near Tier Buffer`,
      size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.nearZeroPadded = new Float32Array(size / 4);

    if (this.midTierBuffer) this.midTierBuffer.destroy();
    this.midTierBuffer = this.device.createBuffer({
      label: `Mid Tier Buffer`,
      size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.midZeroPadded = new Float32Array(size / 4);

    if (this.farTierBuffer) this.farTierBuffer.destroy();
    this.farTierBuffer = this.device.createBuffer({
      label: `Far Tier Buffer`,
      size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.farZeroPadded = new Float32Array(size / 4);

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
      sphereData[f_base + 0] = body.position[0];
      sphereData[f_base + 1] = body.position[1];
      sphereData[f_base + 2] = body.position[2];
      sphereData[f_base + 3] = body.radius;

      // Material
      sphereData.set(body.albedo, f_base + 4); // material.albedo
      // Atmosphere flag in albedo vec4.w
      sphereData[f_base + 7] = body.terrain && (body.terrain as any).atmosphere ? 1.0 : 0.0;
      const emissive = body.emissive ?? [0, 0, 0];
      sphereData.set(emissive, f_base + 8); // material.emissive
      // Planet shading flag: 1.0 if terrain data exists (use detailed terrain shader)
      sphereData[f_base + 11] = body.terrain ? 1.0 : 0.0;
      

      // Planet data
      if (body.terrain) {
        sphereData[f_base + 16] = body.terrain.radius;      // base_radius
        sphereData[f_base + 17] = body.terrain.seaLevel;    // sea_level
        sphereData[f_base + 18] = body.terrain.maxHeight;   // max_height
        sphereData[f_base + 19] = body.terrain.noiseSeed;   // seed
      } else {
        sphereData[f_base + 16] = 0.0;
        sphereData[f_base + 17] = 0.0;
        sphereData[f_base + 18] = 0.0;
        sphereData[f_base + 19] = 0.0;
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


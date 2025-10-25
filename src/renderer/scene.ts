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

    // Write full zero-padded buffers each frame to avoid stale GPU data
    const writeFullBuffer = (buffer: GPUBuffer, bodies: Body[]) => {
      const data = this.serializeSystemState(bodies);
      const zeroPadded = new Float32Array(buffer.size / 4);
      zeroPadded.set(data);
      this.device.queue.writeBuffer(buffer, 0, zeroPadded);
    };

    writeFullBuffer(this.nearTierBuffer, nearBodies);
    writeFullBuffer(this.midTierBuffer, midBodies);
    writeFullBuffer(this.farTierBuffer, farBodies);
  }

  public initializeTierBuffers(initialBodyCount: number) {
    const size = Math.max(1, initialBodyCount * FLOATS_PER_SPHERE * 4);

    if (this.nearTierBuffer) this.nearTierBuffer.destroy();
    this.nearTierBuffer = this.device.createBuffer({
      label: `Near Tier Buffer`,
      size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

    if (this.midTierBuffer) this.midTierBuffer.destroy();
    this.midTierBuffer = this.device.createBuffer({
      label: `Mid Tier Buffer`,
      size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

    if (this.farTierBuffer) this.farTierBuffer.destroy();
    this.farTierBuffer = this.device.createBuffer({
      label: `Far Tier Buffer`,
      size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

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
    const sphereDataU32 = new Uint32Array(sphereData.buffer);
    bodies.forEach((body, i) => {
      const f_base = i * FLOATS_PER_SPHERE;
      const u_base = f_base;
      sphereData[f_base + 0] = body.position[0];
      sphereData[f_base + 1] = body.position[1];
      sphereData[f_base + 2] = body.position[2];
      sphereData[f_base + 3] = body.radius;

      // Material
      sphereData.set(body.albedo, f_base + 4); // material.albedo
      const emissive = body.emissive ?? [0, 0, 0];
      sphereData.set(emissive, f_base + 8); // material.emissive
      sphereData[f_base + 11] = 0.0; // material.fuzziness
      sphereData[f_base + 12] = 1.0; // material.refraction_index

      // Planet data
      if (body.terrain) {
        sphereData[f_base + 16] = body.terrain.radius;      // base_radius
        sphereData[f_base + 17] = body.terrain.seaLevel;    // sea_level
        sphereData[f_base + 18] = body.terrain.maxHeight;   // max_height
        sphereData[f_base + 19] = body.terrain.noiseSeed;   // seed
        sphereDataU32[u_base + 20] = 1; // is_planet = true
      } else {
        sphereData[f_base + 16] = 0.0;
        sphereData[f_base + 17] = 0.0;
        sphereData[f_base + 18] = 0.0;
        sphereData[f_base + 19] = 0.0;
        sphereDataU32[u_base + 20] = 0; // is_planet = false
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


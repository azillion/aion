import type { Body, SystemState, Star } from '../shared/types';
import type { WebGPUCore } from './core';
import type { Camera } from '../camera';

const FLOATS_PER_SPHERE = 16;
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
  public spheresBuffer!: GPUBuffer;
  public starBuffer!: GPUBuffer;

  // A SINGLE, shared buffer for all camera data.
  public sharedCameraUniformBuffer!: GPUBuffer;

  public lastKnownBodyCount: number = 0;
  
  public get hierarchy(): Map<string, string | null> { return this._hierarchy; }
  private _hierarchy: Map<string, string | null> = new Map();

  constructor(core: WebGPUCore) {
    this.device = core.device;
    this.sharedCameraUniformBuffer = this.device.createBuffer({
      label: 'Shared Camera Uniform Buffer',
      size: 256,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }

  public initialize(systemState: SystemState, stars: Star[]) {
    this.lastKnownBodyCount = systemState.bodies.length;
    this.recreateSpheresBuffer(this.lastKnownBodyCount);

    const starData = this.serializeStars(stars);
    this.starBuffer = this.device.createBuffer({
      size: starData.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(this.starBuffer, 0, starData);
  }

  public recreateSpheresBuffer(numBodies: number) {
    if (this.spheresBuffer) {
      this.spheresBuffer.destroy();
    }
    const size = Math.max(1, numBodies * FLOATS_PER_SPHERE * 4);
    this.spheresBuffer = this.device.createBuffer({
      label: `Spheres Buffer (${numBodies} bodies)`,
      size,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    this.lastKnownBodyCount = numBodies;
  }

  public update(systemState: SystemState, camera: Camera, bodiesToRender: Body[], renderScale: number, useWorldSpace: boolean = false) {
    if (systemState.bodies.length !== this.lastKnownBodyCount) {
      this.recreateSpheresBuffer(systemState.bodies.length);
    }
    
    this._hierarchy = buildSystemHierarchy(systemState.bodies);

    const sphereData = this.serializeSystemState(bodiesToRender, camera, renderScale, useWorldSpace);
    this.device.queue.writeBuffer(this.spheresBuffer, 0, sphereData);
  }

  private serializeSystemState(bodies: Body[], camera: Camera, renderScale: number, useWorldSpace: boolean) {
    const sphereData = new Float32Array(this.lastKnownBodyCount * FLOATS_PER_SPHERE);
    const sphereDataU32 = new Uint32Array(sphereData.buffer);
    bodies.forEach((body, i) => {
      const f_base = i * FLOATS_PER_SPHERE;
      const u_base = f_base;
      const renderedRadius = body.radius * renderScale;

      // Sphere.center
      if (useWorldSpace) {
        // For Map and Ship modes, positions are already correctly transformed (scaled or camera-relative)
        sphereData[f_base + 0] = body.position[0];
        sphereData[f_base + 1] = body.position[1];
        sphereData[f_base + 2] = body.position[2];
      } else {
        // Fallback path for modes without preprocessing (kept for future use)
        sphereData[f_base + 0] = body.position[0] - camera.eye[0];
        sphereData[f_base + 1] = body.position[1] - camera.eye[1];
        sphereData[f_base + 2] = body.position[2] - camera.eye[2];
      }
      sphereData[f_base + 3] = renderedRadius; // Sphere.radius

      // Material
      sphereData.set(body.albedo, f_base + 4); // material.albedo
      const emissive = body.emissive ?? [0, 0, 0];
      sphereData.set(emissive, f_base + 8); // material.emissive
      sphereData[f_base + 11] = 0.0; // material.fuzziness
      sphereData[f_base + 12] = 1.0; // material.refraction_index
      sphereDataU32[u_base + 13] = 0; // material.mat_type (Lambertian)
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


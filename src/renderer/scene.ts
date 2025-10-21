import type { Body, SystemState, Vec3, Star } from '../shared/types';
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

  // Camera-related buffers
  public cameraUniformBuffer!: GPUBuffer;
  public galaxyCameraUniformBuffer!: GPUBuffer;

  public lastKnownBodyCount: number = 0;
  
  public get hierarchy(): Map<string, string | null> { return this._hierarchy; }
  private _hierarchy: Map<string, string | null> = new Map();

  constructor(core: WebGPUCore) {
    this.device = core.device;
    this.cameraUniformBuffer = this.device.createBuffer({
      size: 256,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.galaxyCameraUniformBuffer = this.device.createBuffer({
      size: 128,
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

  public update(systemState: SystemState, camera: Camera, bodiesToRender: Body[], renderScale: number) {
    if (systemState.bodies.length !== this.lastKnownBodyCount) {
      this.recreateSpheresBuffer(systemState.bodies.length);
    }
    
    this._hierarchy = buildSystemHierarchy(systemState.bodies);

    const sphereData = this.serializeSystemState(bodiesToRender, camera, renderScale);
    this.device.queue.writeBuffer(this.spheresBuffer, 0, sphereData);

    // Update camera uniform buffer (for compute pass)
    const cameraData = new Float32Array(16);
    cameraData.set([0, 0, 0], 0); // Eye is at origin in camera-relative space
    const lookAtRelative: Vec3 = [
      camera.look_at[0] - camera.eye[0],
      camera.look_at[1] - camera.eye[1],
      camera.look_at[2] - camera.eye[2],
    ];
    cameraData.set(lookAtRelative, 4);
    cameraData.set(camera.up, 8);
    this.device.queue.writeBuffer(this.cameraUniformBuffer, 0, cameraData);
  }

  private serializeSystemState(bodies: Body[], camera: Camera, renderScale: number) {
    const sphereData = new Float32Array(this.lastKnownBodyCount * FLOATS_PER_SPHERE);
    const sphereDataU32 = new Uint32Array(sphereData.buffer);

    bodies.forEach((body, i) => {
      const f_base = i * FLOATS_PER_SPHERE;
      const u_base = f_base;

      const renderedRadius = body.radius * renderScale;

      // Sphere.center (camera-relative)
      sphereData[f_base + 0] = (body.position[0] * renderScale) - camera.eye[0];
      sphereData[f_base + 1] = (body.position[1] * renderScale) - camera.eye[1];
      sphereData[f_base + 2] = (body.position[2] * renderScale) - camera.eye[2];
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
  
  public static calculateRenderScale(bodies: Body[]): number {
    const systemRadius = bodies.reduce((maxDist, b) => {
        const dist = Math.hypot(b.position[0], b.position[1], b.position[2]);
        return Math.max(maxDist, dist);
    }, 0);
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


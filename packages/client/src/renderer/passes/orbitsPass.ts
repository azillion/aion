import { vec3 } from 'gl-matrix';
import type { Body, Theme, Vec3 } from '@shared/types';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import orbitsShaderWGSL from '../shaders/orbits.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import { createShaderModule } from '../shaderUtils';
import { calculateAnalyticOrbit, type OrbitalElements, ORBIT_MAX_POINTS } from '../../physics/orbitalMechanics';


export interface OrbitGlyphData {
  periapsisPoints: (Vec3 | null)[];
  apoapsisPoints: (Vec3 | null)[];
  ascendingNodePoints: (Vec3 | null)[];
  descendingNodePoints: (Vec3 | null)[];
}

const MAX_ORBITS = 256;

export class OrbitsPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private bindGroup!: GPUBindGroup;
  private allOrbitsVertexBuffer!: GPUBuffer; // storage buffer for all orbit points
  private instanceDataBuffer!: GPUBuffer;    // storage buffer for per-instance orbit data
  private visibleInstanceCount = 0;
  private orbitalElements: (OrbitalElements | null)[] = [];
  
  public async initialize(core: WebGPUCore, scene: Scene): Promise<void> {
    const module = await createShaderModule(
        core.device,
        'Orbits Shader Module',
        orbitsShaderWGSL,
        { 'camera.wgsl': cameraWGSL }
    );
    const bindGroupLayout = core.device.createBindGroupLayout({
      label: 'Orbits Bind Group Layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.VERTEX, buffer: { type: 'uniform' } },
        { binding: 1, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } }, // all points
        { binding: 2, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'read-only-storage' } }, // per-instance data
      ],
    });
    this.pipeline = core.device.createRenderPipeline({
      label: 'Orbits Pipeline',
      layout: core.device.createPipelineLayout({
        bindGroupLayouts: [bindGroupLayout],
      }),
      vertex: { module, entryPoint: 'vertexMain' },
      fragment: {
        module, entryPoint: 'fragmentMain',
        targets: [{
          format: core.presentationFormat,
          blend: {
            color: { srcFactor: 'one', dstFactor: 'one', operation: 'add' },
            alpha: { srcFactor: 'one', dstFactor: 'one', operation: 'add' },
          },
        }],
      },
      primitive: { topology: 'triangle-strip' },
    });

    // Allocate large shared buffers
    const pointsPerOrbit = ORBIT_MAX_POINTS;
    const bytesPerPoint = 4 * 4; // vec4<f32>
    const totalPoints = MAX_ORBITS * pointsPerOrbit;
    this.allOrbitsVertexBuffer = core.device.createBuffer({
      size: totalPoints * bytesPerPoint,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

    // Instance data: 8 floats (color vec3 + pointCount + a + e + currentTrueAnomaly) + 2 u32 (offset, reserved)
    const bytesPerInstance = 16 * 4; // pad to 64 bytes for alignment
    this.instanceDataBuffer = core.device.createBuffer({
      size: MAX_ORBITS * bytesPerInstance,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

    this.bindGroup = core.device.createBindGroup({
      label: 'Orbits Bind Group',
      layout: bindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: scene.sharedCameraUniformBuffer } },
        { binding: 1, resource: { buffer: this.allOrbitsVertexBuffer } },
        { binding: 2, resource: { buffer: this.instanceDataBuffer } },
      ],
    });

    this.orbitalElements = Array.from({ length: MAX_ORBITS }, () => null);
  }

  // Deprecated but kept as no-op to preserve call sites if any
  private resizeOrbitBuffers(_device: GPUDevice, _count: number) {}

  public _calculateBarycenter(bodies: Body[]): { position: Vec3, velocity: Vec3, totalMass: number } {
    let totalMass = 0;
    const weightedPos: Vec3 = [0, 0, 0];
    const weightedVel: Vec3 = [0, 0, 0];

    for (const body of bodies) {
        totalMass += body.mass;
        weightedPos[0] += body.position[0] * body.mass;
        weightedPos[1] += body.position[1] * body.mass;
        weightedPos[2] += body.position[2] * body.mass;
        weightedVel[0] += body.velocity[0] * body.mass;
        weightedVel[1] += body.velocity[1] * body.mass;
        weightedVel[2] += body.velocity[2] * body.mass;
    }

    if (totalMass === 0) return { position: [0,0,0], velocity: [0,0,0], totalMass: 0 };

    return {
        position: [weightedPos[0] / totalMass, weightedPos[1] / totalMass, weightedPos[2] / totalMass],
        velocity: [weightedVel[0] / totalMass, weightedVel[1] / totalMass, weightedVel[2] / totalMass],
        totalMass,
    };
  }

  public update(bodies: Body[], _scene: Scene, core: WebGPUCore, systemScale: number): OrbitGlyphData {
    
    const barycenter = this._calculateBarycenter(bodies);
    if (barycenter.totalMass === 0) {
      return { periapsisPoints: [], apoapsisPoints: [], ascendingNodePoints: [], descendingNodePoints: [] };
    }

    const periapsisPoints: (Vec3 | null)[] = Array(bodies.length).fill(null);
    const apoapsisPoints: (Vec3 | null)[] = Array(bodies.length).fill(null);
    const ascendingNodePoints: (Vec3 | null)[] = Array(bodies.length).fill(null);
    const descendingNodePoints: (Vec3 | null)[] = Array(bodies.length).fill(null);

    const pointsPerOrbit = ORBIT_MAX_POINTS;
    const bytesPerPoint = 4 * 4;
    const cpuPoints = new Float32Array(MAX_ORBITS * pointsPerOrbit * 4);
    const instanceStrideFloats = 16; // 64 bytes
    const cpuInstance = new Float32Array(MAX_ORBITS * instanceStrideFloats);

    let visibleIndex = 0;
    for (let i = 0; i < Math.min(bodies.length, MAX_ORBITS); i++) {
      const body = bodies[i];
      const r_vec = vec3.subtract(vec3.create(), body.position, barycenter.position) as Vec3;
      const v_vec = vec3.subtract(vec3.create(), body.velocity, barycenter.velocity) as Vec3;
      const elements = calculateAnalyticOrbit(r_vec, v_vec, barycenter.totalMass, systemScale);
      this.orbitalElements[i] = elements;

      if (elements) {
        // copy points into the large CPU buffer at the correct offset
        const pointBase = visibleIndex * pointsPerOrbit * 4;
        cpuPoints.set(elements.points, pointBase);

        // write instance data
        const instBase = visibleIndex * instanceStrideFloats;
        // color (use theme accent in run, but we can set placeholder; will be overridden if needed)
        // store as neutral white; shaders can colorize
        cpuInstance[instBase + 0] = 1.0; // r
        cpuInstance[instBase + 1] = 1.0; // g
        cpuInstance[instBase + 2] = 1.0; // b
        cpuInstance[instBase + 3] = elements.pointCount;
        cpuInstance[instBase + 4] = elements.a * systemScale;
        cpuInstance[instBase + 5] = elements.e;
        cpuInstance[instBase + 6] = elements.currentTrueAnomaly;
        // pack offset (as float for alignment; shader reads as u32 via bitcast or we place into separate u32 slot)
        cpuInstance[instBase + 7] = pointBase; // offset into points array (in vec4 units)

        periapsisPoints[i] = elements.periapsis;
        apoapsisPoints[i] = elements.apoapsis;
        ascendingNodePoints[i] = elements.ascendingNode;
        descendingNodePoints[i] = elements.descendingNode;
        visibleIndex++;
      }
    }
    // upload to GPU
    core.device.queue.writeBuffer(this.allOrbitsVertexBuffer, 0, cpuPoints);
    core.device.queue.writeBuffer(this.instanceDataBuffer, 0, cpuInstance);
    this.visibleInstanceCount = visibleIndex;
    return { periapsisPoints, apoapsisPoints, ascendingNodePoints, descendingNodePoints };
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, theme: Theme): void {
    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: context.orbitsTexture.createView(),
        loadOp: 'clear', storeOp: 'store', clearValue: { r: 0, g: 0, b: 0, a: 0 },
      }]
    });
    pass.setPipeline(this.pipeline);
    // persistent bind group set once
    pass.setBindGroup(0, this.bindGroup);
    if (this.visibleInstanceCount > 0) {
      pass.draw(ORBIT_MAX_POINTS * 2, this.visibleInstanceCount, 0, 0);
    }
    pass.end();
  }
  
  public clearAll() {
    this.orbitalElements.fill(null);
  }
}

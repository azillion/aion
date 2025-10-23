import { vec3 } from 'gl-matrix';
import type { Body, Theme, Vec3 } from '../../shared/types';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import orbitsShaderWGSL from '../shaders/orbits.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import { createShaderModule } from '../shaderUtils';
import { calculateAnalyticOrbit, type OrbitalElements, ORBIT_MAX_POINTS } from '../../physics/orbitalMechanics';

// ORBIT_MAX_POINTS imported from physics/orbitalMechanics

export class OrbitsPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private uniformBuffer!: GPUBuffer;

  private orbitVertexBuffers: GPUBuffer[] = [];
  private orbitalElements: (OrbitalElements | null)[] = [];

  public periapsisPoints: (Vec3 | null)[] = [];
  public apoapsisPoints: (Vec3 | null)[] = [];
  public ascendingNodePoints: (Vec3 | null)[] = [];
  public descendingNodePoints: (Vec3 | null)[] = [];
  
  public async initialize(core: WebGPUCore, scene: Scene): Promise<void> {
    const module = await createShaderModule(
        core.device,
        'Orbits Shader Module',
        orbitsShaderWGSL,
        { 'camera.wgsl': cameraWGSL }
    );
    this.pipeline = core.device.createRenderPipeline({
      label: 'Orbits Pipeline',
      layout: 'auto',
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

    this.uniformBuffer = core.device.createBuffer({
      size: 48,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    this.resizeOrbitBuffers(core.device, scene.lastKnownBodyCount);
  }

  private resizeOrbitBuffers(device: GPUDevice, count: number) {
    this.orbitVertexBuffers.forEach(b => b.destroy());
    this.orbitVertexBuffers = Array.from({ length: count }, () => device.createBuffer({
        size: ORBIT_MAX_POINTS * 4 * 4,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    }));
    this.orbitalElements = Array.from({ length: count }, () => null);
    this.periapsisPoints = Array.from({ length: count }, () => null);
    this.apoapsisPoints = Array.from({ length: count }, () => null);
    this.ascendingNodePoints = Array.from({ length: count }, () => null);
    this.descendingNodePoints = Array.from({ length: count }, () => null);
  }

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

  public update(bodies: Body[], _scene: Scene, core: WebGPUCore, systemScale: number) {
    if (bodies.length !== this.orbitVertexBuffers.length) {
        this.resizeOrbitBuffers(core.device, bodies.length);
    }
    
    const barycenter = this._calculateBarycenter(bodies);
    if (barycenter.totalMass === 0) return;

    for (let i = 0; i < bodies.length; i++) {
      const body = bodies[i];
      const r_vec = vec3.subtract(vec3.create(), body.position, barycenter.position) as Vec3;
      const v_vec = vec3.subtract(vec3.create(), body.velocity, barycenter.velocity) as Vec3;
      const elements = calculateAnalyticOrbit(r_vec, v_vec, barycenter.totalMass, systemScale);
      this.orbitalElements[i] = elements;

      if (elements) {
        core.device.queue.writeBuffer(this.orbitVertexBuffers[i], 0, elements.points.buffer, 0, elements.points.byteLength);
        this.periapsisPoints[i] = elements.periapsis;
        this.apoapsisPoints[i] = elements.apoapsis;
        this.ascendingNodePoints[i] = elements.ascendingNode;
        this.descendingNodePoints[i] = elements.descendingNode;
      }
    }
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, theme: Theme): void {
    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: context.orbitsTexture.createView(),
        loadOp: 'clear', storeOp: 'store', clearValue: { r: 0, g: 0, b: 0, a: 0 },
      }]
    });
    pass.setPipeline(this.pipeline);

    for (let i = 0; i < this.orbitalElements.length; i++) {
      const elements = this.orbitalElements[i];
      if (!elements || elements.pointCount < 2) continue;

      const uniformData = new Float32Array(12);
      uniformData.set(theme.accent, 0);
      uniformData[3] = elements.pointCount;
      uniformData[4] = elements.a * context.systemScale;
      uniformData[5] = elements.e;
      uniformData[6] = elements.currentTrueAnomaly;
      context.core.device.queue.writeBuffer(this.uniformBuffer, 0, uniformData);

      const bindGroup = context.core.device.createBindGroup({
        label: `Orbits Bind Group ${i}`,
        layout: this.pipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: context.scene.sharedCameraUniformBuffer } },
            { binding: 1, resource: { buffer: this.uniformBuffer } },
            { binding: 2, resource: { buffer: this.orbitVertexBuffers[i] } }
        ],
      });

      pass.setBindGroup(0, bindGroup);
      pass.draw(elements.pointCount * 2);
    }
    pass.end();
  }
  
  public clearAll() {
    this.orbitalElements.fill(null);
  }
}

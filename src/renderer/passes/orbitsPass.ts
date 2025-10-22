import { mat4, vec3 } from 'gl-matrix';
import { G } from '../../shared/constants';
import type { Body, Theme, Vec3 } from '../../shared/types';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import orbitsShaderWGSL from '../shaders/orbits.wgsl?raw';

const ORBIT_SAMPLES = 256;
const ORBIT_MAX_POINTS = ORBIT_SAMPLES + 1;

export class OrbitsPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private uniformBuffer!: GPUBuffer;

  private orbitVertexBuffers: GPUBuffer[] = [];
  private orbitCPUData: Float32Array[] = [];
  private orbitCounts: number[] = [];
  private orbitalElements: ({ a: number, e: number, currentTrueAnomaly: number } | null)[] = [];

  public periapsisPoints: (Vec3 | null)[] = [];
  public apoapsisPoints: (Vec3 | null)[] = [];
  public ascendingNodePoints: (Vec3 | null)[] = [];
  public descendingNodePoints: (Vec3 | null)[] = [];
  

  public initialize(core: WebGPUCore, scene: Scene): void {
    const module = core.device.createShaderModule({ code: orbitsShaderWGSL });
    this.pipeline = core.device.createRenderPipeline({
      label: 'Orbits Pipeline',
      layout: 'auto',
      vertex: {
        module,
        entryPoint: 'vertexMain',
        // No vertex buffers; positions are sourced from a storage buffer in the shader
      },
      fragment: {
        module,
        entryPoint: 'fragmentMain',
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
      size: 112, // WGSL struct min binding size (mat4 + vec3 + f32 + a + e + currentTrueAnomaly + padding)
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    this.resizeOrbitBuffers(core.device, scene.lastKnownBodyCount);
  }

  private resizeOrbitBuffers(device: GPUDevice, count: number) {
    this.orbitVertexBuffers.forEach(b => b.destroy());
    this.orbitVertexBuffers = Array.from({ length: count }, () => device.createBuffer({
        // 16-byte per point (vec4<f32>) to match WGSL storage alignment
        size: ORBIT_MAX_POINTS * 4 * 4,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    }));
    this.orbitCPUData = Array.from({ length: count }, () => new Float32Array(ORBIT_MAX_POINTS * 4));
    this.orbitCounts = Array.from({ length: count }, () => 0);
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
      
      const r_vec: Vec3 = [
        body.position[0] - barycenter.position[0],
        body.position[1] - barycenter.position[1],
        body.position[2] - barycenter.position[2]
      ];
      const v_vec: Vec3 = [
        body.velocity[0] - barycenter.velocity[0],
        body.velocity[1] - barycenter.velocity[1],
        body.velocity[2] - barycenter.velocity[2]
      ];

      this._calculateAnalyticOrbit(r_vec, v_vec, barycenter.totalMass, i, systemScale);
    }
    
    for (let i = 0; i < bodies.length; i++) {
        const bytes = this.orbitCounts[i] * 4 * 4;
        if (bytes > 0) {
            core.device.queue.writeBuffer(this.orbitVertexBuffers[i], 0, this.orbitCPUData[i].buffer, 0, bytes);
        }
    }
  }

  private _calculateAnalyticOrbit(r_vec: Vec3, v_vec: Vec3, totalSystemMass: number, i: number, systemScale: number) {
    const r = Math.hypot(...r_vec);
    const v = Math.hypot(...v_vec);
    const mu = G * totalSystemMass;
    const h_vec = [r_vec[1] * v_vec[2] - r_vec[2] * v_vec[1], r_vec[2] * v_vec[0] - r_vec[0] * v_vec[2], r_vec[0] * v_vec[1] - r_vec[1] * v_vec[0]];
    const h = Math.hypot(...h_vec) || 1e-9;
    const c_vec = [v_vec[1] * h_vec[2] - v_vec[2] * h_vec[1], v_vec[2] * h_vec[0] - v_vec[0] * h_vec[2], v_vec[0] * h_vec[1] - v_vec[1] * h_vec[0]];
    const e_vec = [c_vec[0] / mu - r_vec[0] / r, c_vec[1] / mu - r_vec[1] / r, c_vec[2] / mu - r_vec[2] / r];
    const e = Math.hypot(...e_vec);
    const invA = 2 / r - (v * v) / mu;
    // Robustness: skip non-elliptical or numerically unstable cases (e >= 1)
    if (!isFinite(invA) || invA <= 0 || e >= 1) {
      this.orbitCounts[i] = 0;
      this.orbitalElements[i] = null;
      return;
    }
    const a = 1 / invA;
    const bSemi = a * Math.sqrt(Math.max(0, 1 - e * e));

    const W = [h_vec[0]/h, h_vec[1]/h, h_vec[2]/h];
    let P = e > 1e-6 ? [e_vec[0]/e, e_vec[1]/e, e_vec[2]/e] : [r_vec[0]/r, r_vec[1]/r, r_vec[2]/r];
    const dotPW = P[0]*W[0] + P[1]*W[1] + P[2]*W[2];
    P = [P[0] - dotPW*W[0], P[1] - dotPW*W[1], P[2] - dotPW*W[2]];
    const pLen = Math.hypot(...P) || 1;
    P = [P[0]/pLen, P[1]/pLen, P[2]/pLen];
    const Q = [W[1]*P[2] - W[2]*P[1], W[2]*P[0] - W[0]*P[2], W[0]*P[1] - W[1]*P[0]];

    // Generate orbit ribbon points
    const arr = this.orbitCPUData[i];
    for (let k = 0; k < ORBIT_MAX_POINTS; k++) {
      const theta = (k / (ORBIT_MAX_POINTS - 1)) * Math.PI * 2;
      const xPeri = a * (Math.cos(theta) - e);
      const yPeri = bSemi * Math.sin(theta);
      const base = k * 4;
      arr[base+0] = (P[0]*xPeri + Q[0]*yPeri) * systemScale;
      arr[base+1] = (P[1]*xPeri + Q[1]*yPeri) * systemScale;
      arr[base+2] = (P[2]*xPeri + Q[2]*yPeri) * systemScale;
      arr[base+3] = 1.0; // padding for 16-byte alignment
    }
    this.orbitCounts[i] = ORBIT_MAX_POINTS;

    // Store minimal orbital elements required by shader and glyph pass
    const dotR_Q = r_vec[0]*Q[0] + r_vec[1]*Q[1] + r_vec[2]*Q[2];
    const dotR_P = r_vec[0]*P[0] + r_vec[1]*P[1] + r_vec[2]*P[2];
    const currentNu = Math.atan2(dotR_Q, dotR_P);
    this.orbitalElements[i] = { a, e, currentTrueAnomaly: currentNu };

    // Calculate and store apsis points
    if (e > 1e-6) {
      const Pvec = vec3.fromValues(P[0], P[1], P[2]);
      const periapsisVec = vec3.scale(vec3.create(), Pvec, a * (1.0 - e));
      const apoapsisVec = vec3.scale(vec3.create(), Pvec, -a * (1.0 + e));
      this.periapsisPoints[i] = [
          periapsisVec[0] * systemScale,
          periapsisVec[1] * systemScale,
          periapsisVec[2] * systemScale,
      ];
      this.apoapsisPoints[i] = [
          apoapsisVec[0] * systemScale,
          apoapsisVec[1] * systemScale,
          apoapsisVec[2] * systemScale,
      ];
    } else {
      this.periapsisPoints[i] = null;
      this.apoapsisPoints[i] = null;
    }

    // Calculate and store node points (intersection with ecliptic XY-plane)
    const nodeLine: Vec3 = [-h_vec[1], h_vec[0], 0];
    const nodeLineMag = vec3.length(vec3.fromValues(nodeLine[0], nodeLine[1], nodeLine[2]));
    if (nodeLineMag > 1e-9) {
      const n = vec3.normalize(vec3.create(), vec3.fromValues(nodeLine[0], nodeLine[1], nodeLine[2]));
      const nu = Math.atan2(vec3.dot(n, vec3.fromValues(Q[0], Q[1], Q[2])), vec3.dot(n, vec3.fromValues(P[0], P[1], P[2])));
      const r_dist = a * (1 - e * e) / (1 + e * Math.cos(nu));
      const ascPos = vec3.scale(vec3.create(), n, r_dist * systemScale);
      this.ascendingNodePoints[i] = [ascPos[0], ascPos[1], ascPos[2]];
      const descPos = vec3.scale(vec3.create(), n, -r_dist * systemScale);
      this.descendingNodePoints[i] = [descPos[0], descPos[1], descPos[2]];
    } else {
      this.ascendingNodePoints[i] = null;
      this.descendingNodePoints[i] = null;
    }
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, theme: Theme): void {
    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: context.orbitsTexture.createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: { r: 0, g: 0, b: 0, a: 0 },
      }]
    });
    pass.setPipeline(this.pipeline);
    for (let i = 0; i < this.orbitVertexBuffers.length; i++) {
      if (this.orbitCounts[i] < 2) continue;

      const elements = this.orbitalElements[i];
      if (!elements) continue;

      const uniformData = new Float32Array(28);
      const viewProj = mat4.multiply(mat4.create(), context.camera.projectionMatrix, context.camera.viewMatrix);
      uniformData.set(viewProj, 0);
      uniformData.set(theme.accent, 16);
      // store pointCount at index 19 (float after color vec3)
      uniformData[19] = this.orbitCounts[i];
      uniformData[20] = elements.a * context.systemScale;
      uniformData[21] = elements.e;
      uniformData[22] = elements.currentTrueAnomaly;
      context.core.device.queue.writeBuffer(this.uniformBuffer, 0, uniformData);

      const bindGroup = context.core.device.createBindGroup({
        label: 'Orbits Bind Group',
        layout: this.pipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: this.uniformBuffer } },
            { binding: 1, resource: { buffer: this.orbitVertexBuffers[i] } }
        ],
      });

      pass.setBindGroup(0, bindGroup);
      // Draw two vertices per orbit point; last point will pair with previous via central diff
      pass.draw(this.orbitCounts[i] * 2);
    }
    pass.end();
  }
  
  public clearAll() {
    this.orbitCounts.fill(0);
  }

  
}

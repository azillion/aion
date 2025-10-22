import { mat4, vec3 } from 'gl-matrix';
import { G } from '../../shared/constants';
import type { Body, Orbit, Theme, Vec3 } from '../../shared/types';
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
      },
      fragment: {
        module,
        entryPoint: 'fragmentMain',
        targets: [{
          format: core.presentationFormat,
          blend: {
            color: { srcFactor: 'src-alpha', dstFactor: 'one-minus-src-alpha', operation: 'add' },
            alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add' },
          },
        }],
      },
      primitive: { topology: 'triangle-strip' },
    });

    this.uniformBuffer = core.device.createBuffer({
      size: 112,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    this.resizeOrbitBuffers(core.device, scene.lastKnownBodyCount);
  }

  public resizeOrbitBuffers(device: GPUDevice, count: number) {
    this.orbitVertexBuffers.forEach(b => b.destroy());
    this.orbitVertexBuffers = Array.from({ length: count }, () => device.createBuffer({
        size: ORBIT_MAX_POINTS * 4 * 4,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    }));
    this.periapsisPoints = Array.from({ length: count }, () => null);
    this.apoapsisPoints = Array.from({ length: count }, () => null);
    this.ascendingNodePoints = Array.from({ length: count }, () => null);
    this.descendingNodePoints = Array.from({ length: count }, () => null);
  }

  /**
   * Calculates stable orbital parameters for a body around its parent.
   */
  public static calculateOrbit(body: Body, parent: Body, timestampSeconds: number): Orbit | null {
    const r_vec: Vec3 = [
      body.position[0] - parent.position[0],
      body.position[1] - parent.position[1],
      body.position[2] - parent.position[2],
    ];
    const v_vec: Vec3 = [
      body.velocity[0] - parent.velocity[0],
      body.velocity[1] - parent.velocity[1],
      body.velocity[2] - parent.velocity[2],
    ];
    const r = Math.hypot(...r_vec);
    const v = Math.hypot(...v_vec);
    if (r === 0) return null;

    const mu = G * (parent.mass + body.mass);
    const h_vec = vec3.cross(vec3.create(), r_vec as unknown as vec3, v_vec as unknown as vec3);

    const e_vec_part1 = vec3.cross(vec3.create(), v_vec as unknown as vec3, h_vec);
    const e_vec = vec3.subtract(
      vec3.create(),
      vec3.scale(e_vec_part1, e_vec_part1, 1 / mu),
      vec3.scale(vec3.create(), r_vec as unknown as vec3, 1 / r)
    );
    const eccentricity = vec3.length(e_vec);

    const invA = 2 / r - (v * v) / mu;
    if (!isFinite(invA) || invA <= 0) return null;
    const semiMajorAxis = 1 / invA;

    const W = vec3.normalize(vec3.create(), h_vec);
    let P = eccentricity > 1e-6 ? vec3.normalize(vec3.create(), e_vec) : vec3.normalize(vec3.create(), r_vec as unknown as vec3);
    const dotPW = vec3.dot(P, W);
    P = vec3.normalize(P, vec3.subtract(P, P, vec3.scale(vec3.create(), W, dotPW)));
    const Q = vec3.cross(vec3.create(), W, P);

    const dotR_Q = vec3.dot(r_vec as unknown as vec3, Q);
    const dotR_P = vec3.dot(r_vec as unknown as vec3, P);
    const trueAnomalyAtEpoch = Math.atan2(dotR_Q, dotR_P);

    return {
      semiMajorAxis,
      eccentricity,
      p: [P[0], P[1], P[2]],
      q: [Q[0], Q[1], Q[2]],
      mu,
      trueAnomalyAtEpoch,
      epoch: timestampSeconds,
    };
  }

  /**
   * Returns the position in the parent's frame at the requested absolute time (seconds).
   */
  public static getPositionOnOrbit(orbit: Orbit, timeSeconds: number): Vec3 {
    const dt = timeSeconds - orbit.epoch;
    const n = Math.sqrt(orbit.mu / Math.pow(orbit.semiMajorAxis, 3));
    const M0 = orbit.trueAnomalyAtEpoch - orbit.eccentricity * Math.sin(orbit.trueAnomalyAtEpoch);
    const M = M0 + n * dt;
    let E = M;
    for (let i = 0; i < 5; i++) {
      E = E - (E - orbit.eccentricity * Math.sin(E) - M) / (1 - orbit.eccentricity * Math.cos(E));
    }
    const nu = 2 * Math.atan2(
      Math.sqrt(1 + orbit.eccentricity) * Math.sin(E / 2),
      Math.sqrt(Math.max(0, 1 - orbit.eccentricity)) * Math.cos(E / 2)
    );
    const r_dist = orbit.semiMajorAxis * (1 - orbit.eccentricity * Math.cos(E));

    const x = r_dist * Math.cos(nu);
    const y = r_dist * Math.sin(nu);
    const { p, q } = orbit;
    return [
      p[0] * x + q[0] * y,
      p[1] * x + q[1] * y,
      p[2] * x + q[2] * y,
    ];
  }

  /**
   * Generates or refreshes the orbit ribbon geometry for a body index.
   */
  public updateOrbitGeometry(core: WebGPUCore, bodyIndex: number, orbit: Orbit) {
    const a = orbit.semiMajorAxis;
    const e = orbit.eccentricity;
    const bSemi = a * Math.sqrt(Math.max(0, 1 - e * e));
    const arr = new Float32Array(ORBIT_MAX_POINTS * 4);
    const { p, q } = orbit;
    for (let k = 0; k < ORBIT_MAX_POINTS; k++) {
      const theta = (k / (ORBIT_MAX_POINTS - 1)) * Math.PI * 2;
      const xPeri = a * (Math.cos(theta) - e);
      const yPeri = bSemi * Math.sin(theta);
      const base = k * 4;
      arr[base + 0] = (p[0] * xPeri + q[0] * yPeri);
      arr[base + 1] = (p[1] * xPeri + q[1] * yPeri);
      arr[base + 2] = (p[2] * xPeri + q[2] * yPeri);
      arr[base + 3] = 1.0;
    }
    core.device.queue.writeBuffer(this.orbitVertexBuffers[bodyIndex], 0, arr);

    // Update apsis points for glyphs (nodes optional)
    if (e > 1e-6) {
      const peri = vec3.scale(vec3.create(), vec3.fromValues(orbit.p[0], orbit.p[1], orbit.p[2]), a * (1.0 - e));
      const apo = vec3.scale(vec3.create(), vec3.fromValues(orbit.p[0], orbit.p[1], orbit.p[2]), -a * (1.0 + e));
      this.periapsisPoints[bodyIndex] = [peri[0], peri[1], peri[2]];
      this.apoapsisPoints[bodyIndex] = [apo[0], apo[1], apo[2]];
    } else {
      this.periapsisPoints[bodyIndex] = null;
      this.apoapsisPoints[bodyIndex] = null;
    }
    this.ascendingNodePoints[bodyIndex] = null;
    this.descendingNodePoints[bodyIndex] = null;
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, theme: Theme, orbitsToDraw: (Orbit | null)[], targetTimeSeconds: number): void {
    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: context.orbitsTexture.createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: { r: 0, g: 0, b: 0, a: 0 },
      }]
    });
    pass.setPipeline(this.pipeline);
    // Build a dedicated top-down viewProjection to frame all orbits
    const systemRadius = orbitsToDraw.reduce((max, orbit) => {
      if (!orbit) return max;
      const apo = orbit.semiMajorAxis * (1 + orbit.eccentricity);
      return Math.max(max, apo);
    }, 0);
    if (systemRadius <= 0) { pass.end(); return; }
    const aspect = context.textureSize.width / context.textureSize.height;
    const margin = 1.1;
    const left = -systemRadius * margin * aspect;
    const right = systemRadius * margin * aspect;
    const bottom = -systemRadius * margin;
    const top = systemRadius * margin;
    const near = -systemRadius * 10.0;
    const far = systemRadius * 10.0;
    const proj = mat4.ortho(mat4.create(), left, right, bottom, top, near, far);
    const view = mat4.lookAt(mat4.create(), [0, 0, systemRadius], [0, 0, 0], [0, 1, 0]);
    const viewProj = mat4.multiply(mat4.create(), proj, view);

    for (let i = 0; i < this.orbitVertexBuffers.length; i++) {
      const orbit = orbitsToDraw[i];
      if (!orbit) continue;

      const posOnOrbit = OrbitsPass.getPositionOnOrbit(orbit, targetTimeSeconds);
      const pVec = vec3.fromValues(orbit.p[0], orbit.p[1], orbit.p[2]);
      const qVec = vec3.fromValues(orbit.q[0], orbit.q[1], orbit.q[2]);
      const nu = Math.atan2(vec3.dot(vec3.fromValues(posOnOrbit[0], posOnOrbit[1], posOnOrbit[2]), qVec), vec3.dot(vec3.fromValues(posOnOrbit[0], posOnOrbit[1], posOnOrbit[2]), pVec));

      const uniformData = new Float32Array(28);
      uniformData.set(viewProj, 0);
      uniformData.set(theme.accent, 16);
      uniformData[19] = ORBIT_MAX_POINTS;
      uniformData[20] = orbit.semiMajorAxis;
      uniformData[21] = orbit.eccentricity;
      uniformData[22] = nu;
      context.core.device.queue.writeBuffer(this.uniformBuffer, 0, uniformData);

      const bindGroup = context.core.device.createBindGroup({
        label: 'Orbits Bind Group',
        layout: this.pipeline.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: { buffer: this.uniformBuffer } },
          { binding: 1, resource: { buffer: this.orbitVertexBuffers[i] } },
        ],
      });
      pass.setBindGroup(0, bindGroup);
      pass.draw(ORBIT_MAX_POINTS * 2);
    }
    pass.end();
  }
  
  public clearAll() {
    // No-op with cached approach; geometry will be refreshed as needed
  }
  
  /**
   * Calculates apsides and nodes from a stable Orbit (scaled to systemScale).
   */
  public static getSpecialPoints(orbit: Orbit, systemScale: number): {
    periapsis: Vec3 | null,
    apoapsis: Vec3 | null,
    ascending: Vec3 | null,
    descending: Vec3 | null,
  } {
    if (!isFinite(orbit.semiMajorAxis) || orbit.eccentricity >= 1.0) {
      return { periapsis: null, apoapsis: null, ascending: null, descending: null };
    }
    const pVec = vec3.fromValues(orbit.p[0], orbit.p[1], orbit.p[2]);
    const qVec = vec3.fromValues(orbit.q[0], orbit.q[1], orbit.q[2]);

    const periVec = vec3.scale(vec3.create(), pVec, orbit.semiMajorAxis * (1 - orbit.eccentricity) * systemScale);
    const apoVec = vec3.scale(vec3.create(), pVec, -orbit.semiMajorAxis * (1 + orbit.eccentricity) * systemScale);
    const periapsis: Vec3 = [periVec[0], periVec[1], periVec[2]];
    const apoapsis: Vec3 = [apoVec[0], apoVec[1], apoVec[2]];

    // Node line = intersection of orbital plane (p,q) with ecliptic (Z plane)
    const W = vec3.normalize(vec3.create(), vec3.cross(vec3.create(), pVec, qVec));
    const nodeLine = vec3.fromValues(-W[1], W[0], 0);
    const nodeLen = vec3.length(nodeLine);
    if (nodeLen < 1e-6) {
      return { periapsis, apoapsis, ascending: null, descending: null };
    }
    vec3.scale(nodeLine, nodeLine, 1 / nodeLen);
    const nu = Math.atan2(vec3.dot(nodeLine, qVec), vec3.dot(nodeLine, pVec));
    const r_dist = orbit.semiMajorAxis * (1 - orbit.eccentricity * orbit.eccentricity) / (1 + orbit.eccentricity * Math.cos(nu));
    const asc = vec3.scale(vec3.create(), nodeLine, r_dist * systemScale);
    const dsc = vec3.scale(vec3.create(), nodeLine, -r_dist * systemScale);
    const ascending: Vec3 = [asc[0], asc[1], asc[2]];
    const descending: Vec3 = [dsc[0], dsc[1], dsc[2]];
    return { periapsis, apoapsis, ascending, descending };
  }
}

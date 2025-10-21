import { mat4 } from 'gl-matrix';
import { G } from '../../shared/constants';
import type { Body, Theme } from '../../shared/types';
import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import orbitsShaderWGSL from '../shaders/orbits.wgsl?raw';

const ORBIT_SAMPLES = 256;
const ORBIT_MAX_POINTS = ORBIT_SAMPLES + 1;

export class OrbitsPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private uniformBuffer!: GPUBuffer;
  private bindGroup!: GPUBindGroup;

  private orbitVertexBuffers: GPUBuffer[] = [];
  private orbitCPUData: Float32Array[] = [];
  private orbitCounts: number[] = [];

  public initialize(core: WebGPUCore, scene: Scene): void {
    const module = core.device.createShaderModule({ code: orbitsShaderWGSL });
    this.pipeline = core.device.createRenderPipeline({
      label: 'Orbits Pipeline',
      layout: 'auto',
      vertex: {
        module,
        entryPoint: 'vertexMain',
        buffers: [{
          arrayStride: 12, // 3 * 4 bytes for vec3<f32>
          attributes: [{ shaderLocation: 0, format: 'float32x3', offset: 0 }],
        }],
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
      primitive: { topology: 'line-strip' },
    });

    this.uniformBuffer = core.device.createBuffer({
      size: 80, // mat4x4 + vec3 + padding
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    this.bindGroup = core.device.createBindGroup({
      label: 'Orbits Bind Group',
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [{ binding: 0, resource: { buffer: this.uniformBuffer } }],
    });

    this.resizeOrbitBuffers(core.device, scene.lastKnownBodyCount);
  }

  private resizeOrbitBuffers(device: GPUDevice, count: number) {
    this.orbitVertexBuffers.forEach(b => b.destroy());
    this.orbitVertexBuffers = Array.from({ length: count }, () => device.createBuffer({
        size: ORBIT_MAX_POINTS * 3 * 4,
        usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
    }));
    this.orbitCPUData = Array.from({ length: count }, () => new Float32Array(ORBIT_MAX_POINTS * 3));
    this.orbitCounts = Array.from({ length: count }, () => 0);
  }

  public update(bodies: Body[], scene: Scene, core: WebGPUCore, systemScale: number) {
    if (bodies.length !== this.orbitVertexBuffers.length) {
        this.resizeOrbitBuffers(core.device, bodies.length);
    }
    this.calculateAnalyticOrbits(bodies, scene.hierarchy, systemScale);
    
    for (let i = 0; i < bodies.length; i++) {
        const bytes = this.orbitCounts[i] * 3 * 4;
        if (bytes > 0) {
            core.device.queue.writeBuffer(this.orbitVertexBuffers[i], 0, this.orbitCPUData[i].buffer, 0, bytes);
        }
    }
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, theme: Theme): void {
    const proj = mat4.create();
    const aspect = context.textureSize.height > 0 ? (context.textureSize.width / context.textureSize.height) : 16/9;
    mat4.perspective(proj, Math.PI * 25 / 180, aspect, 0.00001, 10000.0);
    const view = mat4.create();
    mat4.lookAt(view, context.camera.eye as any, context.camera.look_at as any, context.camera.up as any);
    const viewProj = mat4.multiply(mat4.create(), proj, view);
    
    const uniformData = new Float32Array(20);
    uniformData.set(viewProj, 0);
    uniformData.set(theme.accent, 16);
    context.core.device.queue.writeBuffer(this.uniformBuffer, 0, uniformData);

    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: context.orbitsTexture.createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: { r: 0, g: 0, b: 0, a: 0 },
      }]
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, this.bindGroup);
    for (let i = 0; i < this.orbitVertexBuffers.length; i++) {
      if (this.orbitCounts[i] < 2) continue;
      pass.setVertexBuffer(0, this.orbitVertexBuffers[i]);
      pass.draw(this.orbitCounts[i]);
    }
    pass.end();
  }
  
  public clearAll() {
    this.orbitCounts.fill(0);
  }

  private calculateAnalyticOrbits(bodies: Body[], hierarchy: Map<string, string | null>, systemScale: number) {
    for (let i = 0; i < bodies.length; i++) {
      const body = bodies[i];
      const parentId = hierarchy.get(body.id);
      if (!parentId) { this.orbitCounts[i] = 0; continue; }
      const parent = bodies.find(b => b.id === parentId);
      if (!parent) { this.orbitCounts[i] = 0; continue; }

      const r_vec = [body.position[0] - parent.position[0], body.position[1] - parent.position[1], body.position[2] - parent.position[2]];
      const v_vec = [body.velocity[0] - parent.velocity[0], body.velocity[1] - parent.velocity[1], body.velocity[2] - parent.velocity[2]];
      const r = Math.hypot(...r_vec);
      const v = Math.hypot(...v_vec);
      const mu = G * (parent.mass + body.mass);
      const h_vec = [r_vec[1] * v_vec[2] - r_vec[2] * v_vec[1], r_vec[2] * v_vec[0] - r_vec[0] * v_vec[2], r_vec[0] * v_vec[1] - r_vec[1] * v_vec[0]];
      const h = Math.hypot(...h_vec) || 1e-9;
      const c_vec = [v_vec[1] * h_vec[2] - v_vec[2] * h_vec[1], v_vec[2] * h_vec[0] - v_vec[0] * h_vec[2], v_vec[0] * h_vec[1] - v_vec[1] * h_vec[0]];
      const e_vec = [c_vec[0] / mu - r_vec[0] / r, c_vec[1] / mu - r_vec[1] / r, c_vec[2] / mu - r_vec[2] / r];
      const e = Math.hypot(...e_vec);
      const invA = 2 / r - (v * v) / mu;
      if (!isFinite(invA) || invA <= 0) { this.orbitCounts[i] = 0; continue; }
      const a = 1 / invA;
      const bSemi = a * Math.sqrt(Math.max(0, 1 - e * e));
      const W = [h_vec[0]/h, h_vec[1]/h, h_vec[2]/h];
      let P = e > 1e-6 ? [e_vec[0]/e, e_vec[1]/e, e_vec[2]/e] : [r_vec[0]/r, r_vec[1]/r, r_vec[2]/r];
      const dotPW = P[0]*W[0] + P[1]*W[1] + P[2]*W[2];
      P = [P[0] - dotPW*W[0], P[1] - dotPW*W[1], P[2] - dotPW*W[2]];
      const pLen = Math.hypot(...P) || 1;
      P = [P[0]/pLen, P[1]/pLen, P[2]/pLen];
      const Q = [W[1]*P[2] - W[2]*P[1], W[2]*P[0] - W[0]*P[2], W[0]*P[1] - W[1]*P[0]];
      
      const arr = this.orbitCPUData[i];
      for (let k = 0; k < ORBIT_MAX_POINTS; k++) {
        const theta = (k / (ORBIT_MAX_POINTS - 1)) * Math.PI * 2;
        const xPeri = a * (Math.cos(theta) - e);
        const yPeri = bSemi * Math.sin(theta);
        const base = k * 3;
        arr[base+0] = (parent.position[0] + P[0]*xPeri + Q[0]*yPeri) * systemScale;
        arr[base+1] = (parent.position[1] + P[1]*xPeri + Q[1]*yPeri) * systemScale;
        arr[base+2] = (parent.position[2] + P[2]*xPeri + Q[2]*yPeri) * systemScale;
      }
      this.orbitCounts[i] = ORBIT_MAX_POINTS;
    }
  }
}

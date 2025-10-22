import type { Body, Vec3 } from './shared/types';
import type { Camera } from './camera';
import { CameraMode } from './state';
import { mat4, vec4 } from 'gl-matrix';

export class HUDManager {
  private canvas: HTMLCanvasElement;
  private context: CanvasRenderingContext2D;

  constructor(canvas: HTMLCanvasElement) {
    this.canvas = canvas;
    const ctx = this.canvas.getContext('2d');
    if (!ctx) throw new Error('Failed to get 2D context for HUD');
    this.context = ctx;
  }
  
  public getCanvas(): HTMLCanvasElement {
    return this.canvas;
  }

  public draw(bodies: Body[], camera: Camera, viewport: { width: number, height: number }, cameraMode: CameraMode, playerShipId: string | null): void {
    this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);

    // Style
    this.context.font = '12px Inter, system-ui, Avenir, Helvetica, Arial, sans-serif';
    this.context.textAlign = 'left';
    this.context.textBaseline = 'top';
    this.context.fillStyle = 'rgba(255,255,255,0.9)';
    this.context.strokeStyle = 'rgba(255,255,255,0.9)';

    if (cameraMode === CameraMode.GALACTIC_MAP) {
      return;
    }

    const centerX = viewport.width * 0.5;
    const centerY = viewport.height * 0.5;

    let iterable: Body[] = bodies;
    if (cameraMode === CameraMode.SHIP_RELATIVE && playerShipId) {
      iterable = bodies.filter(b => b.id !== playerShipId);
    }

    for (const body of iterable) {
      const p = this._projectToScreen(body.position, camera, viewport);
      if (!p) continue;

      const isFocus = camera.focusBodyId === body.id;
      if (cameraMode === CameraMode.SYSTEM_MAP) {
        if (!p.inView) continue;
        const x = p.x;
        const y = p.y;
        this.context.beginPath();
        this.context.arc(x, y, isFocus ? 4 : 2.5, 0, Math.PI * 2);
        this.context.fill();
        this.context.fillText(body.name, x + 6, y - 6);
        continue;
      }

      // SHIP_RELATIVE: show on-screen and off-screen indicators (no labels)
      if (p.inView) {
        const x = p.x;
        const y = p.y;
        this.context.beginPath();
        this.context.arc(x, y, isFocus ? 4 : 3, 0, Math.PI * 2);
        this.context.fill();
      } else {
        const x = p.x;
        const y = p.y;
        // Draw a chevron pointing towards the off-screen target
        const angle = Math.atan2(y - centerY, x - centerX);
        const size = isFocus ? 10 : 8;
        this._drawChevron(x, y, angle, size);
      }
    }
  }

  private _projectToScreen(pos: Vec3, camera: Camera, viewport: { width: number, height: number }): { x: number, y: number, inView: boolean } | null {
    if (camera.isOrthographic) {
      const viewProj = mat4.multiply(mat4.create(), camera.projectionMatrix as unknown as number[], camera.viewMatrix as unknown as number[]);
      const worldPos = vec4.fromValues(pos[0], pos[1], pos[2], 1.0);
      const clip = vec4.transformMat4(vec4.create(), worldPos, viewProj as unknown as number[]);
      const w = clip[3];
      if (!Number.isFinite(w) || w === 0) return null;
      const ndcX = clip[0] / w;
      const ndcY = clip[1] / w;
      let px = (ndcX * 0.5 + 0.5) * viewport.width;
      let py = (-ndcY * 0.5 + 0.5) * viewport.height;
      return { x: px, y: py, inView: true };
    }
    const eye = camera.eye as unknown as number[];
    const look = camera.look_at as unknown as number[];
    const dx = eye[0] - look[0];
    const dy = eye[1] - look[1];
    const dz = eye[2] - look[2];
    const dist = Math.max(1e-9, Math.hypot(dx, dy, dz));

    const vfov = 25 * Math.PI / 180;
    const aspect = (viewport.height > 0) ? (viewport.width / viewport.height) : (16 / 9);
    const near = Math.max(1e-6, dist * 0.001);
    const far = 1e15;
    const proj = mat4.create();
    mat4.perspective(proj, vfov, aspect, near, far);

    const view = mat4.create();
    mat4.lookAt(view, eye, look, camera.up as unknown as number[]);

    const viewProj = mat4.multiply(mat4.create(), proj, view);

    const worldPos = vec4.fromValues(pos[0], pos[1], pos[2], 1.0);
    const clip = vec4.transformMat4(vec4.create(), worldPos, viewProj);
    const w = clip[3];
    if (!Number.isFinite(w) || w <= 0) return null;
    const ndcX = clip[0] / w;
    const ndcY = clip[1] / w;
    let px = (ndcX * 0.5 + 0.5) * viewport.width;
    let py = (-ndcY * 0.5 + 0.5) * viewport.height;

    const inView = px >= 0 && px <= viewport.width && py >= 0 && py <= viewport.height;
    if (inView) return { x: px, y: py, inView: true };

    // Clamp to screen edge along ray from center to (px, py)
    const cx = viewport.width * 0.5;
    const cy = viewport.height * 0.5;
    const vx = px - cx;
    const vy = py - cy;
    const eps = 1e-9;
    let tx = Number.POSITIVE_INFINITY;
    let ty = Number.POSITIVE_INFINITY;
    if (Math.abs(vx) > eps) {
      const tx1 = (0 - cx) / vx;
      const tx2 = ((viewport.width - 1) - cx) / vx;
      tx = vx > 0 ? tx2 : tx1;
    }
    if (Math.abs(vy) > eps) {
      const ty1 = (0 - cy) / vy;
      const ty2 = ((viewport.height - 1) - cy) / vy;
      ty = vy > 0 ? ty2 : ty1;
    }
    const t = Math.min(tx, ty);
    px = cx + vx * t;
    py = cy + vy * t;
    px = Math.max(0, Math.min(viewport.width - 1, px));
    py = Math.max(0, Math.min(viewport.height - 1, py));
    return { x: px, y: py, inView: false };
  }

  private _drawChevron(x: number, y: number, angle: number, size: number): void {
    this.context.save();
    this.context.translate(x, y);
    this.context.rotate(angle);
    this.context.beginPath();
    this.context.moveTo(size, 0);
    this.context.lineTo(-size * 0.6, size * 0.5);
    this.context.lineTo(-size * 0.6, -size * 0.5);
    this.context.closePath();
    this.context.fill();
    this.context.restore();
  }
}

import type { Camera } from '../camera/camera';
import type { Vec3 } from '../shared/types';
import { vec4 } from 'gl-matrix';

export function projectWorldToScreen(
  worldPos: Vec3,
  camera: Camera,
  viewport: { width: number; height: number }
): { x: number; y: number; inView: boolean } | null {
  const viewProj = camera.viewProjectionMatrix;
  const posVec4 = vec4.fromValues(worldPos[0], worldPos[1], worldPos[2], 1.0);

  const clip = vec4.transformMat4(vec4.create(), posVec4, viewProj as unknown as number[]);

  const w = clip[3];
  if (!Number.isFinite(w) || w <= 0) {
    // Point is behind the camera's near plane or at infinity.
    return null;
  }

  const ndcX = clip[0] / w;
  const ndcY = clip[1] / w;

  // Check if in view using NDC (-1 to +1 range is inside the view frustum)
  const inView = ndcX >= -1 && ndcX <= 1 && ndcY >= -1 && ndcY <= 1;

  // Convert from Normalized Device Coordinates (-1 to +1) to screen coordinates (0 to width/height).
  let px = (ndcX * 0.5 + 0.5) * viewport.width;
  let py = (-ndcY * 0.5 + 0.5) * viewport.height; // Y is inverted in NDC

  if (inView) {
    return { x: px, y: py, inView: true };
  }

  // If not in view, we still return a coordinate clamped to the screen edge.
  // This is useful for off-screen indicators.
  const cx = viewport.width * 0.5;
  const cy = viewport.height * 0.5;
  const vx = px - cx;
  const vy = py - cy;
  const eps = 1e-9;
  let tx = Number.POSITIVE_INFINITY;
  let ty = Number.POSITIVE_INFINITY;

  if (Math.abs(vx) > eps) {
    const tx1 = (0 - cx) / vx;
    const tx2 = (viewport.width - 1 - cx) / vx;
    tx = vx > 0 ? tx2 : tx1;
  }
  if (Math.abs(vy) > eps) {
    const ty1 = (0 - cy) / vy;
    const ty2 = (viewport.height - 1 - cy) / vy;
    ty = vy > 0 ? ty2 : ty1;
  }

  const t = Math.min(tx, ty);
  px = cx + vx * t;
  py = cy + vy * t;

  // Final clamp to ensure it's on the border.
  px = Math.max(0, Math.min(viewport.width - 1, px));
  py = Math.max(0, Math.min(viewport.height - 1, py));

  return { x: px, y: py, inView: false };
}



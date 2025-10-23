import type { Body, FrameData, Vec3 } from './shared/types';
import { CameraMode } from './state';
import { projectWorldToScreen } from './renderer/projection';

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

  public draw(frameData: FrameData): void {
    const { bodiesToRender, camera, viewport, cameraMode, playerShipId } = frameData;
    this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);

    // Style
    this.context.font = '12px Inter, system-ui, Avenir, Helvetica, Arial, sans-serif';
    this.context.textAlign = 'left';
    this.context.textBaseline = 'top';
    this.context.fillStyle = 'rgba(255,255,255,0.9)';
    this.context.strokeStyle = 'rgba(255,255,255,0.9)';

    // Do not draw any HUD elements in galactic map or ship view
    if (cameraMode === CameraMode.GALACTIC_MAP || cameraMode === CameraMode.SHIP_RELATIVE) {
      return;
    }

    const centerX = viewport.width * 0.5;
    const centerY = viewport.height * 0.5;

    const iterable: Body[] = bodiesToRender;

    for (const body of iterable) {
      const p = projectWorldToScreen(body.position as Vec3, camera, viewport);
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

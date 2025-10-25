import type { FrameData, Vec3 } from './shared/types';
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

  public clear(): void {
    this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);
  }

  public draw(frameData: FrameData): void {
    const { bodiesToRender, camera, viewport, cameraMode } = frameData;
    this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);

    // Style
    this.context.font = '12px Inter, system-ui, Avenir, Helvetica, Arial, sans-serif';
    this.context.textAlign = 'left';
    this.context.textBaseline = 'top';
    this.context.fillStyle = 'rgba(255,255,255,0.9)';
    this.context.strokeStyle = 'rgba(255,255,255,0.9)';

    // Do not draw any HUD elements in galactic map
    if (cameraMode === CameraMode.GALACTIC_MAP) {
      return;
    }

    // const centerX = viewport.width * 0.5;
    // const centerY = viewport.height * 0.5;

    // --- New Ship Relative HUD Logic ---
    if (cameraMode === CameraMode.SHIP_RELATIVE) {
      this.drawShipHUD(frameData);
      // We return here because drawShipHUD handles its own body indicators
      return;
    }
    // --- End New Logic ---

    for (const body of bodiesToRender) {
      const p = projectWorldToScreen(body.position as Vec3, camera, viewport);
      if (!p) continue;

      const isFocus = camera.focusBodyId === body.id;
      // This block is now only for SYSTEM_MAP
      if (cameraMode === CameraMode.SYSTEM_MAP && p.inView) {
        this.context.beginPath();
        this.context.arc(p.x, p.y, isFocus ? 4 : 2.5, 0, Math.PI * 2);
        this.context.fill();
        this.context.fillText(body.name, p.x + 6, p.y - 6);
      }
    }
  }

  

  // --- Ship Relative HUD ---
  private drawShipHUD(frameData: FrameData): void {
    const { rawState, camera, viewport, playerShipId } = frameData;
    const playerShip = rawState.bodies.find(b => b.id === playerShipId);
    if (!playerShip || !('velocity' in playerShip)) {
      return;
    }
    const ship = playerShip as import('./shared/types').Ship;

    const centerX = viewport.width / 2;
    const centerY = viewport.height / 2;

    this.context.strokeStyle = 'rgba(148, 255, 179, 0.7)';
    this.context.fillStyle = 'rgba(148, 255, 179, 0.7)';
    this.context.lineWidth = 1.5;

    // 1. Draw Boresight (where the nose is pointing)
    this.drawBoresight(centerX, centerY);

    // 2. Draw Velocity Vector (Prograde Marker)
    // For simplicity, we'll calculate prograde relative to the sun (solar system barycenter)
    const sun = rawState.bodies.find(b => b.id === 'sol');
    const sunVel = sun ? sun.velocity : [0,0,0];
    const relativeVel: Vec3 = [
      ship.velocity[0] - sunVel[0],
      ship.velocity[1] - sunVel[1],
      ship.velocity[2] - sunVel[2]
    ];
    const speed = Math.hypot(...relativeVel);

    // Only draw the prograde marker if we are moving
    const showPrograde = false;
    if (showPrograde && speed > 1) { // 1 km/s threshold
      // The prograde vector is just the normalized velocity.
      // We can scale it by an arbitrary amount to make it visible.
      // This position is ALREADY relative to the ship (which is at the origin).
      const progradeRelativePos: Vec3 = [
        relativeVel[0] / speed * 1000, // Scaled for visibility
        relativeVel[1] / speed * 1000,
        relativeVel[2] / speed * 1000
      ];
      const screenPos = projectWorldToScreen(progradeRelativePos, camera, viewport);
      if (screenPos && screenPos.inView) {
        this.drawProgradeMarker(screenPos.x, screenPos.y, false);
      }
    }

    // 3. Draw Speed Readout
    this.context.font = '16px "Share Tech Mono", monospace';
    this.context.textAlign = 'center';
    this.context.fillText(`${(speed).toFixed(1)} km/s`, centerX, centerY + 40);

    // 4. Target Info
    const focusBody = rawState.bodies.find(b => b.id === camera.focusBodyId);
    if (focusBody && focusBody.id !== playerShipId) {
      const dist = Math.hypot(
        focusBody.position[0] - ship.position[0],
        focusBody.position[1] - ship.position[1],
        focusBody.position[2] - ship.position[2]
      );
      this.context.textAlign = 'right';
      this.context.fillText(`TARGET: ${focusBody.name.toUpperCase()}`, viewport.width - 12, viewport.height - 30);
      this.context.fillText(`RANGE: ${(dist/1_000_000).toFixed(2)} Gm`, viewport.width - 12, viewport.height - 12);
    }

    // 5. Indicators: Precision / Kill Rotation
    const flags = rawState.flags || {};
    // Place indicators at edges: top-right for Precision, top-left for Kill Rot
    if (flags.precision) {
      this.context.textAlign = 'right';
      this.context.fillText('PRECISION', viewport.width - 12, 12);
    }
    if (flags.killRotation) {
      this.context.textAlign = 'left';
      this.context.fillText('KILL ROT', 12, 12);
    }

    // 6. Body Indicators: Draw indicators for off-screen bodies
    this.drawBodyIndicators(frameData);
  }

  private drawBoresight(x: number, y: number): void {
    const size = 10;
    this.context.beginPath();
    this.context.moveTo(x + size, y);
    this.context.lineTo(x + size / 2, y);
    this.context.moveTo(x - size, y);
    this.context.lineTo(x - size / 2, y);
    this.context.moveTo(x, y + size);
    this.context.lineTo(x, y + size / 2);
    this.context.moveTo(x, y - size);
    this.context.lineTo(x, y - size / 2);
    this.context.stroke();
  }

  private drawProgradeMarker(x: number, y: number, isRetrograde: boolean): void {
    const radius = 12;
    this.context.beginPath();
    this.context.arc(x, y, radius, 0, Math.PI * 2);
    this.context.moveTo(x - radius, y);
    this.context.lineTo(x - radius - 5, y);
    this.context.moveTo(x + radius, y);
    this.context.lineTo(x + radius + 5, y);
    this.context.moveTo(x, y - radius);
    this.context.lineTo(x, y - radius - 5);
    
    if (!isRetrograde) {
      this.context.moveTo(x, y + radius);
      this.context.lineTo(x, y + radius + 5);
    }
    this.context.stroke();
  }

  private drawBodyIndicators(frameData: FrameData) {
    // This is the logic from the original 'draw' method, now separated
    const { bodiesToRender, camera, viewport, playerShipId } = frameData;
    const centerX = viewport.width * 0.5;
    const centerY = viewport.height * 0.5;
    
    const iterable = bodiesToRender.filter(b => b.id !== playerShipId);

    for (const body of iterable) {
      const p = projectWorldToScreen(body.position as Vec3, camera, viewport);
      if (!p) continue;
      
      const isFocus = camera.focusBodyId === body.id;
      if (p.inView) {
        this.context.beginPath();
        this.context.arc(p.x, p.y, isFocus ? 4 : 3, 0, Math.PI * 2);
        this.context.fill();
      } else {
        const angle = Math.atan2(p.y - centerY, p.x - centerX);
        const size = isFocus ? 10 : 8;
        this._drawChevron(p.x, p.y, angle, size);
      }
    }
  }
  // --- End Ship Relative HUD ---
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

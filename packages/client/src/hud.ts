import type { Vec3 } from '@shared/types';
import { CameraMode } from './state';
import { projectWorldToScreen } from './renderer/projection';
import type { RenderPayload } from './views/types';

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
    this.context.save();
    this.context.setTransform(1, 0, 0, 1, 0, 0);
    this.context.clearRect(0, 0, this.canvas.width, this.canvas.height);
    this.context.restore();
  }

  public draw(payload: RenderPayload): void {
    this.clear();

    // Style
    this.context.font = '12px "Share Tech Mono", monospace';
    this.context.textAlign = 'left';
    this.context.textBaseline = 'top';
    this.context.fillStyle = 'rgba(255,255,255,0.9)';
    this.context.strokeStyle = 'rgba(255,255,255,0.9)';
    this.context.lineWidth = 1.0;
    this.context.lineCap = 'butt';
    this.context.lineJoin = 'miter';
    
    if (payload.cameraMode === CameraMode.SHIP_RELATIVE) {
      this.drawShipHUD(payload);
      // We return here because drawShipHUD handles its own body indicators
      return;
    }
  }

  

  // --- Ship Relative HUD ---
  private drawShipHUD(payload: RenderPayload): void {
    const { rawState, camera, viewport, playerShipId } = payload;
    const playerShip = rawState.ships.find(s => s.body.id === playerShipId);
    if (!playerShip || !('velocity' in playerShip.body)) {
      return;
    }
    const ship = playerShip;

    const centerX = viewport.width / 2;
    const centerY = viewport.height / 2;

    this.context.strokeStyle = 'rgba(255, 255, 255, 0.85)';
    this.context.fillStyle = 'rgba(255, 255, 255, 0.85)';
    this.context.lineWidth = 1.0;

    // 1. Draw Boresight (where the nose is pointing)
    this.drawBoresight(centerX, centerY);

    // 2. Draw Velocity Vector (Prograde Marker)
    // For simplicity, we'll calculate prograde relative to the sun (solar system barycenter)
    const sun = rawState.bodies.find(b => b.id === 'sol');
    const sunVel = sun ? sun.velocity : [0,0,0];
    const relativeVel: Vec3 = [
      ship.body.velocity[0] - sunVel[0],
      ship.body.velocity[1] - sunVel[1],
      ship.body.velocity[2] - sunVel[2]
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
    const focusBody = rawState.bodies.find(b => b.id === camera.focusBodyId) ?? rawState.ships.find(s => s.body.id === camera.focusBodyId)?.body;
    if (focusBody && focusBody.id !== playerShipId) {
      const dx = focusBody.position[0] - ship.body.position[0];
      const dy = focusBody.position[1] - ship.body.position[1];
      const dz = focusBody.position[2] - ship.body.position[2];
      const distKm = Math.hypot(dx, dy, dz);

      // Terrain-aware surface radius (km)
      const baseR = (focusBody as any).terrain?.radius ?? focusBody.radius;
      const sea = (focusBody as any).terrain?.seaLevel ?? 0.0;
      const maxH = ((focusBody as any).terrain?.maxHeight ?? 0.0) * baseR;
      const surfaceR = baseR + Math.max(sea, maxH);
      const altitudeKm = Math.max(0, distKm - surfaceR);

      const rangeStr = this._formatDistance(distKm);
      const altStr = this._formatDistance(altitudeKm);

      this.context.textAlign = 'right';
      this.context.fillText(`TARGET: ${focusBody.name.toUpperCase()}`, viewport.width - 12, viewport.height - 44);
      this.context.fillText(`ALT: ${altStr}`, viewport.width - 12, viewport.height - 28);
      this.context.fillText(`RANGE: ${rangeStr}`, viewport.width - 12, viewport.height - 12);
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
    this.drawTargetIndicators(payload);
  }

  public resize(width: number, height: number, dpr: number = (globalThis as any).devicePixelRatio || 1): void {
    // Set the canvas's internal pixel buffer to device pixels for crisp lines
    this.canvas.width = Math.max(1, Math.floor(width * dpr));
    this.canvas.height = Math.max(1, Math.floor(height * dpr));
    // Ensure CSS size matches logical size (so drawing API uses CSS pixels)
    this.canvas.style.width = `${width}px`;
    this.canvas.style.height = `${height}px`;
    // Scale drawing operations so coordinates remain in CSS pixels
    this.context.setTransform(dpr, 0, 0, dpr, 0, 0);
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

  private drawTargetIndicators(payload: RenderPayload) {
    const { bodiesToRender, camera, viewport, playerShipId, worldCameraEye } = payload;
    const centerX = viewport.width * 0.5;
    const centerY = viewport.height * 0.5;
    
    const iterable = bodiesToRender ?? [];

    for (const body of iterable.filter(b => b.id !== playerShipId)) {
      const relativePos: Vec3 = [
        (body.position as Vec3)[0] - worldCameraEye[0],
        (body.position as Vec3)[1] - worldCameraEye[1],
        (body.position as Vec3)[2] - worldCameraEye[2],
      ];
      const p = projectWorldToScreen(relativePos, camera, viewport);
      if (!p) continue;
      
      const isFocus = camera.focusBodyId === body.id;
      if (p.inView) {
        const size = isFocus ? 5 : 3;
        this.context.beginPath();
        this.context.rect(p.x - size, p.y - size, size * 2, size * 2);
        this.context.fill();
        if (isFocus) {
          this.context.strokeRect(p.x - size - 2, p.y - size - 2, (size + 2) * 2, (size + 2) * 2);
        }
      } else {
        const angle = Math.atan2(p.y - centerY, p.x - centerX);
        if (isFocus) {
          const size = 10;
          this._drawChevron(p.x, p.y, angle, size);
        } else {
          // Draw a small dot for other off-screen targets
          const size = 2;
          this.context.beginPath();
          this.context.arc(p.x, p.y, size, 0, Math.PI * 2);
          this.context.fill();
        }
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
    this.context.stroke();
    this.context.restore();
  }

  // Distance formatting helper (input: km)
  private _formatDistance(km: number): string {
    const AU_KM = 149_597_870.7;
    if (km >= 0.1 * AU_KM) {
      const au = km / AU_KM;
      return `${au.toFixed(3)} AU`;
    }
    if (km >= 1_000_000) {
      const gm = km / 1_000_000; // 1 Gm = 1e6 km
      return `${gm.toFixed(gm >= 100 ? 0 : 2)} Gm`;
    }
    if (km >= 1_000) {
      const mm = km / 1_000; // 1 Mm = 1e3 km
      return `${mm.toFixed(mm >= 100 ? 0 : 1)} Mm`;
    }
    if (km >= 1) {
      return `${km.toFixed(km >= 100 ? 0 : km >= 10 ? 1 : 2)} km`;
    }
    const m = km * 1000;
    return `${m.toFixed(m >= 100 ? 0 : m >= 10 ? 1 : 2)} m`;
  }
}

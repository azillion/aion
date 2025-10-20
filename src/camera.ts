import type { Vec3, Body, Ship } from './shared/types';
import { vec3, quat } from 'gl-matrix';

export class Camera {
	public eye: Vec3 = [0, 0.5, 5.0];
	public look_at: Vec3 = [0, 0, 0];
	public up: Vec3 = [0, 1, 0];

	public focusBodyId: string | null = null;
	public pendingFrame = false;

	private canvas: HTMLCanvasElement;
	private isDragging = false;
	private lastMouseX = 0;
	private lastMouseY = 0;
	private minDistance = 1e-6;

	constructor(canvas: HTMLCanvasElement) {
		this.canvas = canvas;
		this.canvas.addEventListener('mousedown', this.onMouseDown);
		this.canvas.addEventListener('mouseup', this.onMouseUp);
		this.canvas.addEventListener('mousemove', this.onMouseMove);
		this.canvas.addEventListener('wheel', this.onWheel, { passive: false } as any);
		window.addEventListener('keydown', this.onKeyDown);
	}

	public setFocus(bodyId: string): void {
		this.focusBodyId = bodyId;
		this.pendingFrame = true;
	}

	public update(bodies: Body[], scale: number): void {
		if (this.focusBodyId === null) return;

		const focusBody = bodies.find(b => b.id === this.focusBodyId);
		if (focusBody) {
			const prevTarget = this.look_at;

			const nextTargetX = focusBody.position[0] * scale;
			const nextTargetY = focusBody.position[1] * scale;
			const nextTargetZ = focusBody.position[2] * scale;

			// Calculate how much the target has moved
			const dx = nextTargetX - prevTarget[0];
			const dy = nextTargetY - prevTarget[1];
			const dz = nextTargetZ - prevTarget[2];

			// Update the look_at point to the new position
			this.look_at[0] = nextTargetX;
			this.look_at[1] = nextTargetY;
			this.look_at[2] = nextTargetZ;

			// Simply translate the camera's eye by the same amount
			this.eye[0] += dx;
			this.eye[1] += dy;
			this.eye[2] += dz;

			// Framing now occurs via completePendingFrame(radiusWorld) invoked by renderer
		}
	}

	public completePendingFrame(radiusWorld: number, opts?: { renderScale?: number, desiredPixelRadius?: number, viewportHeight?: number, vfovDeg?: number }): void {
		if (!this.pendingFrame) return;

		let desiredDist: number;
		if (opts && opts.desiredPixelRadius && opts.viewportHeight && opts.vfovDeg && opts.renderScale) {
			const rRender = Math.max(1e-12, radiusWorld * opts.renderScale);
			const theta = (opts.vfovDeg * Math.PI) / 180.0;
			const projScale = 1.0 / Math.tan(theta / 2.0);
			const pixelsPerNdc = 0.5 * opts.viewportHeight;
			// p_px = (r/z) * projScale * pixelsPerNdc ⇒ z = r * projScale * pixelsPerNdc / p_px
			desiredDist = (rRender * projScale * pixelsPerNdc) / Math.max(1.0, opts.desiredPixelRadius);
		} else {
			const size = (Number.isFinite(radiusWorld) && radiusWorld > 0) ? radiusWorld : 1.0;
			desiredDist = size * 4.0;
		}
		desiredDist = Math.max(this.minDistance, desiredDist);

		// Move along current view direction to preserve orientation
		let dirX = this.eye[0] - this.look_at[0];
		let dirY = this.eye[1] - this.look_at[1];
		let dirZ = this.eye[2] - this.look_at[2];
		let len = Math.hypot(dirX, dirY, dirZ);
		if (len <= 1e-9) { dirX = 0.0; dirY = 0.0; dirZ = 1.0; len = 1.0; }
		const nx = dirX / len;
		const ny = dirY / len;
		const nz = dirZ / len;

		this.eye = [
			this.look_at[0] + nx * desiredDist,
			this.look_at[1] + ny * desiredDist,
			this.look_at[2] + nz * desiredDist,
		];

		this.pendingFrame = false;
	}

	private onMouseDown = (event: MouseEvent) => {
		this.isDragging = true;
		this.lastMouseX = event.clientX;
		this.lastMouseY = event.clientY;
	};

	private onWheel = (event: WheelEvent) => {
		// Exponential zoom across large scales; prevent page scroll
		event.preventDefault();
		const sensitivity = 0.0015; // tune for responsiveness
		const multiplier = Math.exp(event.deltaY * sensitivity);
		this._zoomBy(multiplier);
	};

	private onKeyDown = (event: KeyboardEvent) => {
		if (event.ctrlKey || event.metaKey || event.altKey) return;
		const target = event.target as HTMLElement | null;
		if (target && (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA' || (target as HTMLElement).isContentEditable)) return;
		const key = event.key;
		const code = event.code;
		if (key === '+' || key === '=' || code === 'NumpadAdd') {
			this.zoomIn();
		} else if (key === '-' || key === '_' || code === 'NumpadSubtract') {
			this.zoomOut();
		}
	};

	private onMouseUp = () => {
		this.isDragging = false;
	};

	private onMouseMove = (event: MouseEvent) => {
		if (!this.isDragging) return;

		const dx = event.clientX - this.lastMouseX;
		const dy = event.clientY - this.lastMouseY;

		// Simple yaw/pitch around target without external deps
		const toEyeX = this.eye[0] - this.look_at[0];
		const toEyeY = this.eye[1] - this.look_at[1];
		const toEyeZ = this.eye[2] - this.look_at[2];

		const yaw = -dx * 0.005;
		const pitch = -dy * 0.005;

		// Apply yaw around world up (0,1,0)
		const cosYaw = Math.cos(yaw);
		const sinYaw = Math.sin(yaw);
		let rx = cosYaw * toEyeX + sinYaw * toEyeZ;
		let rz = -sinYaw * toEyeX + cosYaw * toEyeZ;

		// Apply pitch around camera right vector approximated via cross(up, forward)
		const forwardX = -rx;
		const forwardY = -toEyeY;
		const forwardZ = -rz;
		const rightX = this.up[1] * forwardZ - this.up[2] * forwardY;
		const rightY = this.up[2] * forwardX - this.up[0] * forwardZ;
		const rightZ = this.up[0] * forwardY - this.up[1] * forwardX;
		const rightLen = Math.hypot(rightX, rightY, rightZ) || 1;
		const rnx = rightX / rightLen;
		const rny = rightY / rightLen;
		const rnz = rightZ / rightLen;

		const cosPitch = Math.cos(pitch);
		const sinPitch = Math.sin(pitch);
		const dotRP = rnx * rx + rny * toEyeY + rnz * rz;
		rx = rx * cosPitch + (rnx * dotRP) * (1 - cosPitch) + (rny * rz - rnz * toEyeY) * sinPitch;
		const ryTemp = toEyeY * cosPitch + (rny * dotRP) * (1 - cosPitch) + (rnz * rx - rnx * rz) * sinPitch;
		rz = rz * cosPitch + (rnz * dotRP) * (1 - cosPitch) + (rnx * toEyeY - rny * rx) * sinPitch;

		this.eye = [this.look_at[0] + rx, this.look_at[1] + ryTemp, this.look_at[2] + rz];

		this.lastMouseX = event.clientX;
		this.lastMouseY = event.clientY;
	};

	public zoomIn(): void {
		this._zoomBy(0.9);
	}

	public zoomOut(): void {
		this._zoomBy(1.1);
	}

	private _zoomBy(multiplier: number): void {
		const dirX = this.eye[0] - this.look_at[0];
		const dirY = this.eye[1] - this.look_at[1];
		const dirZ = this.eye[2] - this.look_at[2];
		const len = Math.hypot(dirX, dirY, dirZ) || 1;
		const nx = dirX / len;
		const ny = dirY / len;
		const nz = dirZ / len;

		const nextDist = Math.max(this.minDistance, len * multiplier);
		this.eye = [
			this.look_at[0] + nx * nextDist,
			this.look_at[1] + ny * nextDist,
			this.look_at[2] + nz * nextDist,
		];
	}

	public updateShipRelative(ship: Ship): void {
		// First-person: eye at ship position; forward from orientation quaternion
		this.eye = ship.position;
		const q = quat.fromValues(ship.orientation[0], ship.orientation[1], ship.orientation[2], ship.orientation[3]);
		const baseForward = vec3.fromValues(0, 0, -1);
		const baseUp = vec3.fromValues(0, 1, 0);
		const forward = vec3.create();
		const up = vec3.create();
		vec3.transformQuat(forward, baseForward, q);
		vec3.transformQuat(up, baseUp, q);
		this.look_at = [
			ship.position[0] + forward[0],
			ship.position[1] + forward[1],
			ship.position[2] + forward[2],
		];
		this.up = [up[0], up[1], up[2]];
	}
}



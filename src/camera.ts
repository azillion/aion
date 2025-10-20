import type { Vec3, Body } from './shared/types';

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
	private minDistance = 0.005;

	constructor(canvas: HTMLCanvasElement) {
		this.canvas = canvas;
		this.canvas.addEventListener('mousedown', this.onMouseDown);
		this.canvas.addEventListener('mouseup', this.onMouseUp);
		this.canvas.addEventListener('mousemove', this.onMouseMove);
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

	public completePendingFrame(radiusWorld: number): void {
		if (!this.pendingFrame) return;

		let desiredDist = radiusWorld * 4.0;
		desiredDist = Math.max(this.minDistance, desiredDist);

		const viewDirX = 0.5;
		const viewDirY = 0.5;
		const viewDirZ = 1.0;
		const len = Math.hypot(viewDirX, viewDirY, viewDirZ) || 1;

		const nx = viewDirX / len;
		const ny = viewDirY / len;
		const nz = viewDirZ / len;

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
}



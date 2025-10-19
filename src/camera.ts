import type { Vec3, Body } from './shared/types';

export class Camera {
	public eye: Vec3 = [0, 0.3, 1.0];
	public look_at: Vec3 = [0, 0, -2.0];
	public up: Vec3 = [0, 1, 0];

	public focusBodyId: string | null = null;
	private pendingFrame = false;

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
			const prevTargetX = this.look_at[0];
			const prevTargetY = this.look_at[1];
			const prevTargetZ = this.look_at[2];

			const nextTargetX = focusBody.position[0] * scale;
			const nextTargetY = focusBody.position[1] * scale;
			const nextTargetZ = focusBody.position[2] * scale - 2.0;

			const dx = nextTargetX - prevTargetX;
			const dy = nextTargetY - prevTargetY;
			const dz = nextTargetZ - prevTargetZ;

			// Translate camera to preserve relative offset while tracking the moving target
			this.look_at[0] = nextTargetX;
			this.look_at[1] = nextTargetY;
			this.look_at[2] = nextTargetZ;
			this.eye[0] += dx;
			this.eye[1] += dy;
			this.eye[2] += dz;

			// If focus just changed, adjust zoom to fit the focused body on screen
			if (this.pendingFrame) {
				// Match shader's camera vertical FOV (degrees)
				const vfovDeg = 25.0;
				const vfovRad = (vfovDeg * Math.PI) / 180.0;
				const tanHalfV = Math.tan(vfovRad / 2.0) || 1e-6;
				const aspect = this.canvas.height > 0 ? (this.canvas.width / this.canvas.height) : 16 / 9;
				// Match render radius scaling used in _serializeSystemState
				const radiusWorld = (focusBody.radius * scale * 100);
				const padding = 1.2; // slightly back off to avoid clipping
				const distVert = (radiusWorld * padding) / tanHalfV;
				const distHorz = (radiusWorld * padding) / (Math.max(1e-6, aspect) * tanHalfV);
				let desiredDist = Math.max(distVert, distHorz);
				desiredDist = Math.max(this.minDistance, desiredDist);

				const dirX = this.eye[0] - this.look_at[0];
				const dirY = this.eye[1] - this.look_at[1];
				const dirZ = this.eye[2] - this.look_at[2];
				const len = Math.hypot(dirX, dirY, dirZ);
				const nx = (len > 1e-6) ? (dirX / len) : 0;
				const ny = (len > 1e-6) ? (dirY / len) : 0;
				const nz = (len > 1e-6) ? (dirZ / len) : 1; // default to +Z if degenerate

				this.eye = [
					this.look_at[0] + nx * desiredDist,
					this.look_at[1] + ny * desiredDist,
					this.look_at[2] + nz * desiredDist,
				];

				this.pendingFrame = false;
			}
		}
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



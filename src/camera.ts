import type { Vec3 } from './shared/types';

export class Camera {
	public eye: Vec3 = [0, 0.3, 1.0];
	public look_at: Vec3 = [0, 0, -2.0];
	public up: Vec3 = [0, 1, 0];

	private canvas: HTMLCanvasElement;
	private isDragging = false;
	private lastMouseX = 0;
	private lastMouseY = 0;
	private zoom = 1.0;

	constructor(canvas: HTMLCanvasElement) {
		this.canvas = canvas;
		this.canvas.addEventListener('mousedown', this.onMouseDown);
		this.canvas.addEventListener('mouseup', this.onMouseUp);
		this.canvas.addEventListener('mousemove', this.onMouseMove);
		this.canvas.addEventListener('wheel', this.onWheel, { passive: false });
	}

	private onMouseDown = (event: MouseEvent) => {
		this.isDragging = true;
		this.lastMouseX = event.clientX;
		this.lastMouseY = event.clientY;
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

	private onWheel = (event: WheelEvent) => {
		event.preventDefault();
		this.zoom += event.deltaY * -0.001;
		this.zoom = Math.max(0.1, Math.min(this.zoom, 5.0));

		const dirX = this.eye[0] - this.look_at[0];
		const dirY = this.eye[1] - this.look_at[1];
		const dirZ = this.eye[2] - this.look_at[2];
		const len = Math.hypot(dirX, dirY, dirZ) || 1;
		const nx = dirX / len;
		const ny = dirY / len;
		const nz = dirZ / len;

		const dist = len * this.zoom;
		this.eye = [this.look_at[0] + nx * dist, this.look_at[1] + ny * dist, this.look_at[2] + nz * dist];
	};
}



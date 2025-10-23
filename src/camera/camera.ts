import type { Vec3 } from '../shared/types';
import { mat4 } from 'gl-matrix';

export class Camera {
	public eye: Vec3 = [0, 0.5, 5.0];
	public look_at: Vec3 = [0, 0, 0];
	public up: Vec3 = [0, 1, 0];
	public right: Vec3 = [1, 0, 0];
	public forward: Vec3 = [0, 0, -1];

	public focusBodyId: string | null = null;
	public pendingFrame = false;

	public projectionMatrix: mat4 = mat4.create();
	public viewMatrix: mat4 = mat4.create();
	public viewProjectionMatrix: mat4 = mat4.create();
	public isOrthographic: boolean = false;
	public distanceToTarget: number = 0;

	public recomputeViewFromLook(): void {
		mat4.lookAt(this.viewMatrix, this.eye as unknown as number[], this.look_at as unknown as number[], this.up as unknown as number[]);
	}

	public refreshDerived(): void {
		const m = this.viewMatrix as unknown as number[];
		const r: Vec3 = [m[0], m[4], m[8]];
		const u: Vec3 = [m[1], m[5], m[9]];
		const f: Vec3 = [-m[2], -m[6], -m[10]];

		const normalize = (v: Vec3): Vec3 => {
			const len = Math.hypot(v[0], v[1], v[2]) || 1.0;
			return [v[0] / len, v[1] / len, v[2] / len];
		};

		this.right = normalize(r);
		this.up = normalize(u);
		this.forward = normalize(f);

		mat4.multiply(this.viewProjectionMatrix, this.projectionMatrix, this.viewMatrix);

		const dx = this.look_at[0] - this.eye[0];
		const dy = this.look_at[1] - this.eye[1];
		const dz = this.look_at[2] - this.eye[2];
		this.distanceToTarget = Math.hypot(dx, dy, dz);
	}

	public updateViewMatrix(): void {
		this.recomputeViewFromLook();
		this.refreshDerived();
	}

	public clone(): Camera {
		const newCam = new Camera();
		newCam.eye = [this.eye[0], this.eye[1], this.eye[2]];
		newCam.look_at = [this.look_at[0], this.look_at[1], this.look_at[2]];
		newCam.up = [this.up[0], this.up[1], this.up[2]];
		newCam.right = [this.right[0], this.right[1], this.right[2]];
		newCam.forward = [this.forward[0], this.forward[1], this.forward[2]];
		newCam.focusBodyId = this.focusBodyId;
		newCam.isOrthographic = this.isOrthographic;
		newCam.distanceToTarget = this.distanceToTarget;
		mat4.copy(newCam.projectionMatrix, this.projectionMatrix);
		return newCam;
	}
}



import type { Vec3 } from '../shared/types';
import { mat4 } from 'gl-matrix';

export class Camera {
	public eye: Vec3 = [0, 0.5, 5.0];
	public look_at: Vec3 = [0, 0, 0];
	public up: Vec3 = [0, 1, 0];

	public focusBodyId: string | null = null;
	public pendingFrame = false;

	public projectionMatrix: mat4 = mat4.create();
	public viewMatrix: mat4 = mat4.create();
	public isOrthographic: boolean = false;
}



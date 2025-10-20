import type { Vec3 } from '../shared/types';

export class Camera {
	public eye: Vec3 = [0, 0.5, 5.0];
	public look_at: Vec3 = [0, 0, 0];
	public up: Vec3 = [0, 1, 0];

	public focusBodyId: string | null = null;
	public pendingFrame = false;
}



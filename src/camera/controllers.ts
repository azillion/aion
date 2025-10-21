import type { Ship } from '../shared/types';
import type { Camera } from '../camera';
import { vec3, quat, mat4 } from 'gl-matrix';
import type { Body } from '../shared/types';

export interface ICameraController {
	update(camera: Camera, context: any): void;
}

export class SystemMapController implements ICameraController {
    public update(camera: Camera, context: { bodies: Body[], scale: number, viewport: { width: number, height: number } }): void {
        camera.isOrthographic = true;

        const focusBody = context.bodies.find(b => b.id === camera.focusBodyId) || context.bodies[0];
        if (!focusBody) return;

        const center_x = focusBody.position[0] * context.scale;
        const center_y = focusBody.position[1] * context.scale;

        const viewDistance = 30.0;
        const aspect = context.viewport.width / Math.max(1, context.viewport.height);

        mat4.ortho(camera.projectionMatrix, -viewDistance * aspect, viewDistance * aspect, -viewDistance, viewDistance, -100.0, 100.0);
        
        camera.eye = [center_x, center_y, 50];
        camera.look_at = [center_x, center_y, 0];
        camera.up = [0, 1, 0];
        mat4.lookAt(camera.viewMatrix, camera.eye as unknown as number[], camera.look_at as unknown as number[], camera.up as unknown as number[]);
    }
}

export class ShipRelativeController implements ICameraController {
	public update(camera: Camera, context: { playerShip: Ship | undefined | null }): void {
		const ship = context.playerShip;
		if (!ship) return;
		camera.eye = ship.position;
		const q = quat.fromValues(ship.orientation[0], ship.orientation[1], ship.orientation[2], ship.orientation[3]);
		const baseForward = vec3.fromValues(0, 0, -1);
		const baseUp = vec3.fromValues(0, 1, 0);
		const forward = vec3.create();
		const up = vec3.create();
		vec3.transformQuat(forward, baseForward, q);
		vec3.transformQuat(up, baseUp, q);
		camera.look_at = [
			ship.position[0] + forward[0],
			ship.position[1] + forward[1],
			ship.position[2] + forward[2],
		];
		camera.up = [up[0], up[1], up[2]];
	}
}



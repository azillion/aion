import type { Ship } from '../shared/types';
import type { Camera } from '../camera';
import { vec3, quat, mat4 } from 'gl-matrix';
import type { Body } from '../shared/types';
import { ReferenceFrame } from '../state';

export interface ICameraController {
	update(camera: Camera, context: any): void;
}

export class SystemMapController implements ICameraController {
    public update(camera: Camera, context: { bodies: Body[], scale: number, viewport: { width: number, height: number }, referenceFrame: ReferenceFrame }): void {
        camera.isOrthographic = true;

        let center_x = 0;
        let center_y = 0;

        const viewDistance = 30.0;
        const aspect = context.viewport.width / Math.max(1, context.viewport.height);

        mat4.ortho(camera.projectionMatrix, -viewDistance * aspect, viewDistance * aspect, -viewDistance, viewDistance, -100.0, 100.0);
        
        camera.eye = [center_x, center_y, 50];
        camera.look_at = [center_x, center_y, 0];
        camera.up = [0, 1, 0];
        mat4.lookAt(camera.viewMatrix, camera.eye as unknown as number[], camera.look_at as unknown as number[], camera.up as unknown as number[]);
        mat4.rotateX(camera.viewMatrix, camera.viewMatrix, -0.3);
    }
}

export class ShipRelativeController implements ICameraController {
    public update(camera: Camera, context: { playerShip: Ship | undefined | null, targetBody?: Body }): void {
        const ship = context.playerShip;
        if (!ship) return;
        camera.eye = ship.position;

        const q = quat.fromValues(ship.orientation[0], ship.orientation[1], ship.orientation[2], ship.orientation[3]);
        const shipUp = vec3.fromValues(0, 1, 0);
        vec3.transformQuat(shipUp, shipUp, q);

        if (context.targetBody) {
            camera.look_at = [
                context.targetBody.position[0],
                context.targetBody.position[1],
                context.targetBody.position[2],
            ];
        } else {
            const shipForward = vec3.fromValues(0, 0, -1);
            vec3.transformQuat(shipForward, shipForward, q);
            camera.look_at = [
                ship.position[0] + shipForward[0],
                ship.position[1] + shipForward[1],
                ship.position[2] + shipForward[2],
            ];
        }

        // Build a stable orthonormal basis
        const forward = vec3.create();
        vec3.subtract(forward, camera.look_at as unknown as vec3, camera.eye as unknown as vec3);
        vec3.normalize(forward, forward);

        let right = vec3.create();
        vec3.cross(right, forward, shipUp);
        // Stage 1 fallback: use world up if shipUp is collinear with forward
        if (vec3.length(right) < 1e-3) {
            const worldUp = vec3.fromValues(0, 1, 0);
            vec3.cross(right, forward, worldUp);
        }
        // Stage 2 fallback: use world Z as guaranteed non-collinear axis
        if (vec3.length(right) < 1e-3) {
            const worldZ = vec3.fromValues(0, 0, 1);
            vec3.cross(right, forward, worldZ);
        }
        vec3.normalize(right, right);

        const finalUp = vec3.create();
        vec3.cross(finalUp, right, forward);

        camera.up = [finalUp[0], finalUp[1], finalUp[2]];
    }
}



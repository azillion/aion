import type { Ship, Body } from '../shared/types';
import type { Camera } from '../camera';
import { vec3, quat, mat4 } from 'gl-matrix';
import { ReferenceFrame } from '../state';

export interface ICameraController {
	update(camera: Camera, context: any): void;
}

export interface SystemMapControllerContext {
  bodies: Body[];
  scale: number;
  viewport: { width: number; height: number };
  referenceFrame: ReferenceFrame;
}

export interface ShipRelativeControllerContext {
  playerShip: Ship | undefined | null;
  targetBody?: Body;
}

export class SystemMapController implements ICameraController {
    public update(camera: Camera, context: SystemMapControllerContext): void {
        camera.isOrthographic = true;

        let center_x = 0;
        let center_y = 0;

        const viewDistance = 30.0;
        const aspect = context.viewport.width / Math.max(1, context.viewport.height);

        mat4.ortho(camera.projectionMatrix, -viewDistance * aspect, viewDistance * aspect, -viewDistance, viewDistance, -100.0, 100.0);

        camera.eye = [center_x, center_y, 50];
        camera.look_at = [center_x, center_y, 0];
        camera.up = [0, 1, 0];
        camera.recomputeViewFromLook();
        mat4.rotateX(camera.viewMatrix, camera.viewMatrix, -0.3);
        camera.refreshDerived();
    }
}

export class ShipRelativeController implements ICameraController {
    public update(camera: Camera, context: ShipRelativeControllerContext): void {
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

		// The mat4.lookAt function is robust. We only need to provide the ship's
		// local "up" vector as a hint. It will handle cases where the 'forward'
		// and 'up' vectors are nearly collinear and build a stable basis.
		camera.up = [shipUp[0], shipUp[1], shipUp[2]];
        camera.recomputeViewFromLook();
        camera.refreshDerived();
    }
}



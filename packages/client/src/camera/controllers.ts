import type { Ship, Body, Vec3 } from '@shared/types';
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
  keys: Set<string>;
  viewport: { width: number; height: number };
}

export interface ShipRelativeLookAtContext {
  keys: Set<string>;
  relativeBodies: Body[];
  playerShip?: Ship | null;
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
        camera.isOrthographic = false;
        const aspect = context.viewport.width / Math.max(1, context.viewport.height);
        mat4.perspective(camera.projectionMatrix, camera.vfov * (Math.PI / 180.0), aspect, 0.001, 1e10);
        const ship = context.playerShip;
        if (!ship) return;
        camera.eye = ship.body.position;

        const shipOrientation = quat.fromValues(
            ship.orientation[0],
            ship.orientation[1],
            ship.orientation[2],
            ship.orientation[3]
        );

        // Set UP from ship orientation; LOOK_AT handled separately using camera-relative data
        const shipUp = vec3.fromValues(0, 1, 0);
        vec3.transformQuat(shipUp, shipUp, shipOrientation);
        camera.up = [shipUp[0], shipUp[1], shipUp[2]];
    }

    public updateLookAt(camera: Camera, context: ShipRelativeLookAtContext): void {
		const isLookingAtTarget = context.keys.has('KeyT');
		const focusBody = camera.focusBodyId ? context.relativeBodies.find(b => b.id === camera.focusBodyId) : undefined;

		if (isLookingAtTarget && focusBody) {
			// Focused body position is already camera-relative; safe and precise
			camera.look_at = focusBody.position as Vec3;
			return;
		}

        // Default: look forward from the cockpit (camera eye will be at origin by app loop)
        const ship = context.playerShip ?? null;
        const shipForward = vec3.fromValues(0, 0, -1);
        if (ship) {
            const shipOrientation = quat.fromValues(
                ship.orientation[0],
                ship.orientation[1],
                ship.orientation[2],
                ship.orientation[3]
            );
            vec3.transformQuat(shipForward, shipForward, shipOrientation);
        }
        camera.look_at = [shipForward[0], shipForward[1], shipForward[2]];
    }
}



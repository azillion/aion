import type { Body, Ship } from '../shared/types';
import type { Camera } from '../camera';
import { vec3, quat } from 'gl-matrix';

export interface ICameraController {
	update(camera: Camera, context: any): void;
}

export class SystemOrbitalController implements ICameraController {
	private minDistance = 1e-6;
	private cameraRef: Camera | null = null;

	constructor(_canvas: HTMLCanvasElement) {
		window.addEventListener('keydown', this.onKeyDown);
	}

	public update(camera: Camera, context: { bodies: Body[], scale: number, viewport: { width: number, height: number }, vfov: number }): void {
		this.cameraRef = camera;
		const bodies = context.bodies;
		const scale = context.scale;
		const viewportHeight = Math.max(1, context.viewport?.height ?? 1);
		const vfovDeg = context.vfov ?? 25.0;

		if (camera.focusBodyId === null) return;
		const focusBody = bodies.find(b => b.id === camera.focusBodyId);
		if (!focusBody) return;

		const prevTarget = camera.look_at;
		const nextTargetX = focusBody.position[0] * scale;
		const nextTargetY = focusBody.position[1] * scale;
		const nextTargetZ = focusBody.position[2] * scale;

		const dx = nextTargetX - prevTarget[0];
		const dy = nextTargetY - prevTarget[1];
		const dz = nextTargetZ - prevTarget[2];

		camera.look_at[0] = nextTargetX;
		camera.look_at[1] = nextTargetY;
		camera.look_at[2] = nextTargetZ;

		camera.eye[0] += dx;
		camera.eye[1] += dy;
		camera.eye[2] += dz;

		if (camera.pendingFrame) {
			const desiredPixels = 100.0;
			const vfovRad = (vfovDeg * Math.PI) / 180.0;
			const renderedRadius = Math.max(1e-12, focusBody.radius * scale);
			const projScale = 1.0 / Math.tan(vfovRad / 2.0);
			// distance = (r / tan(fov/2)) * (viewport_height / desired_px)
			let desiredDist = (renderedRadius * projScale * (viewportHeight * 0.5)) / Math.max(1.0, desiredPixels * 0.5);
			desiredDist = Math.max(this.minDistance, desiredDist);

			let dirX = camera.eye[0] - camera.look_at[0];
			let dirY = camera.eye[1] - camera.look_at[1];
			let dirZ = camera.eye[2] - camera.look_at[2];
			let len = Math.hypot(dirX, dirY, dirZ);
			if (len <= 1e-9) { dirX = 0.0; dirY = 0.0; dirZ = 1.0; len = 1.0; }
			const nx = dirX / len;
			const ny = dirY / len;
			const nz = dirZ / len;
			camera.eye = [
				camera.look_at[0] + nx * desiredDist,
				camera.look_at[1] + ny * desiredDist,
				camera.look_at[2] + nz * desiredDist,
			];
			camera.pendingFrame = false;
		}
	}


	private onKeyDown = (event: KeyboardEvent) => {
		const camera = this.cameraRef;
		if (!camera) return;
		if (event.ctrlKey || event.metaKey || event.altKey) return;
		const target = event.target as HTMLElement | null;
		if (target && (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA' || (target as HTMLElement).isContentEditable)) return;
		const key = event.key;
		const code = event.code;
		if (key === '+' || key === '=' || code === 'NumpadAdd') {
			this._zoomBy(camera, 0.9);
		} else if (key === '-' || key === '_' || code === 'NumpadSubtract') {
			this._zoomBy(camera, 1.1);
		}
	};

	private _zoomBy(camera: Camera, multiplier: number): void {
		const dirX = camera.eye[0] - camera.look_at[0];
		const dirY = camera.eye[1] - camera.look_at[1];
		const dirZ = camera.eye[2] - camera.look_at[2];
		const len = Math.hypot(dirX, dirY, dirZ) || 1;
		const nx = dirX / len;
		const ny = dirY / len;
		const nz = dirZ / len;
		const minD = 1e-5;
		const maxD = 1e12;
		const unclamped = len * multiplier;
		const nextDist = Math.min(maxD, Math.max(minD, unclamped));
		camera.eye = [
			camera.look_at[0] + nx * nextDist,
			camera.look_at[1] + ny * nextDist,
			camera.look_at[2] + nz * nextDist,
		];
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



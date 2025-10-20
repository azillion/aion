import { CameraMode } from '../state';
import { Camera } from './camera';
import { SystemOrbitalController, ShipRelativeController } from './controllers';
import type { Body, Ship } from '../shared/types';

export class CameraManager {
	private camera: Camera;
	private systemOrbitalController: SystemOrbitalController;
	private shipRelativeController: ShipRelativeController;

	constructor(canvas: HTMLCanvasElement) {
		this.camera = new Camera();
		this.systemOrbitalController = new SystemOrbitalController(canvas);
		this.shipRelativeController = new ShipRelativeController();
	}

	public update(mode: CameraMode, context: any): void {
		switch (mode) {
			case CameraMode.SYSTEM_ORBITAL: {
				this.systemOrbitalController.update(this.camera, context as { bodies: Body[], scale: number, focusBodyRenderedRadius: number });
				break;
			}
			case CameraMode.SHIP_RELATIVE: {
				this.shipRelativeController.update(this.camera, context as { playerShip: Ship | undefined | null });
				break;
			}
			case CameraMode.GALACTIC_MAP: {
				// No-op for now; galactic map uses renderer-built matrices
				break;
			}
			default: break;
		}
	}

	public getCamera(): Camera {
		return this.camera;
	}
}



import { CameraMode } from '../state';
import { Camera } from './camera';
import { SystemMapController, ShipRelativeController } from './controllers';
import type { Ship } from '../shared/types';

export class CameraManager {
	private camera: Camera;
	private systemMapController: SystemMapController;
	private shipRelativeController: ShipRelativeController;

	constructor(_canvas: HTMLCanvasElement) {
		this.camera = new Camera();
		this.systemMapController = new SystemMapController();
		this.shipRelativeController = new ShipRelativeController();
	}

	public update(mode: CameraMode, context: any): void {
		switch (mode) {
			case CameraMode.SYSTEM_MAP: {
				this.systemMapController.update(this.camera, context as any);
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



import { CameraMode } from '../state';
import { Camera } from './camera';
import { SystemMapController, ShipRelativeController, type SystemMapControllerContext, type ShipRelativeControllerContext } from './controllers';

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
                if ('bodies' in context && 'scale' in context) {
                    this.systemMapController.update(this.camera, context as SystemMapControllerContext);
                }
                break;
            }
            case CameraMode.SHIP_RELATIVE: {
                if ('playerShip' in context) {
                    this.shipRelativeController.update(this.camera, context as ShipRelativeControllerContext);
                }
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



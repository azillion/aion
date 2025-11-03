import { CameraMode } from '../state';
import { Camera } from './camera';
import { SystemMapController, ShipRelativeController, type SystemMapControllerContext, type ShipRelativeControllerContext, type ShipRelativeLookAtContext } from './controllers';

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
            case CameraMode.SHIP_RELATIVE: {
                if ('playerShip' in context) {
                    this.shipRelativeController.update(this.camera, context as ShipRelativeControllerContext);
                }
                break;
            }
            default: break;
        }
    }

    public updateLookAt(camera: Camera, context: ShipRelativeLookAtContext): void {
        this.shipRelativeController.updateLookAt(camera, context);
    }

	public getCamera(): Camera {
		return this.camera;
	}
}



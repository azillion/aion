export enum CameraMode {
  SHIP_RELATIVE,
  SYSTEM_MAP,
  GALACTIC_MAP,
}

export class AppState {
  public cameraMode: CameraMode = CameraMode.SHIP_RELATIVE;
  public showOrbits: boolean = false;
  public playerShipId: string | null = null;
  public crtIntensity: number = 1.0;
}



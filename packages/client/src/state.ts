export enum CameraMode {
  SHIP_RELATIVE,
}

export class AppState {
  public cameraMode: CameraMode = CameraMode.SHIP_RELATIVE; // Default and only mode for now
  public showOrbits: boolean = false;
  public showAtmosphere: boolean = true;
  public showHUD: boolean = true;
  public playerShipId: string | null = null;
  public crtIntensity: number = 1.0;
}



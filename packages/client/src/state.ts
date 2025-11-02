export enum ReferenceFrame {
  BARYCENTRIC = 'Barycentric',
  FOCUSED_BODY = 'Focused Body',
}

export enum CameraMode {
  SHIP_RELATIVE,
  SYSTEM_MAP,
}

export class AppState {
  public cameraMode: CameraMode = CameraMode.SHIP_RELATIVE; // Default and only mode for now
  public showOrbits: boolean = false;
  public showAtmosphere: boolean = true;
  public showHUD: boolean = true;
  public playerShipId: string | null = null;
  public crtIntensity: number = 1.0;
  public referenceFrame: ReferenceFrame = ReferenceFrame.BARYCENTRIC;
  public referenceBodyId: string | null = 'sol';
}



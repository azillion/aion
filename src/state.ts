export enum ReferenceFrame {
  BARYCENTRIC = 'Barycentric',
  FOCUSED_BODY = 'Focused Body',
}

export enum CameraMode {
  SHIP_RELATIVE,
  SYSTEM_MAP,
  GALACTIC_MAP,
}

export class AppState {
  public cameraMode: CameraMode = CameraMode.SHIP_RELATIVE;
  public showOrbits: boolean = false;
  public showHUD: boolean = true;
  public playerShipId: string | null = null;
  public crtIntensity: number = 1.0;
  public referenceFrame: ReferenceFrame = ReferenceFrame.BARYCENTRIC;
  public referenceBodyId: string | null = 'sol';
  public debugTierView: number = -1; // -1: off, 0: Near, 1: Mid, 2: Far
}



export enum ViewMode {
  System,
  Galaxy,
}

export class AppState {
  public viewMode: ViewMode = ViewMode.System;
  public showOrbits: boolean = false;
}



declare module 'stats.js' {
  export default class Stats {
    dom: HTMLDivElement;
    showPanel(panel: number): void;
    begin(): void;
    end(): number;
    update(): void;
  }
}



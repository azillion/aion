declare module 'dat.gui' {
  export class GUI {
    add(...args: any[]): any;
    addFolder(...args: any[]): any;
    remove(controller: any): void;
    open(): void;
    close(): void;
    dom?: HTMLElement;
  }
}

declare namespace dat {
  class GUIController {
    updateDisplay(): void;
  }
}



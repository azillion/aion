export class InputManager {
    public deltaX = 0;
    public deltaY = 0;
    public keys: Set<string> = new Set();

    constructor(_canvas: HTMLCanvasElement) {
        document.addEventListener('keydown', (e: KeyboardEvent) => {
            this.keys.add(e.code);
        });
        document.addEventListener('keyup', (e: KeyboardEvent) => {
            this.keys.delete(e.code);
        });
    }

    public tick(): void {
        this.deltaX = 0;
        this.deltaY = 0;
    }
}



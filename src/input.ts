export class InputManager {
    public deltaX = 0;
    public deltaY = 0;
    public keys: Set<string> = new Set();
    private canvas: HTMLCanvasElement;

    constructor(canvas: HTMLCanvasElement) {
        this.canvas = canvas;
        document.addEventListener('mousemove', (e: MouseEvent) => {
            if (document.pointerLockElement === this.canvas) {
                this.deltaX += e.movementX;
                this.deltaY += e.movementY;
            }
        });
        document.addEventListener('keydown', (e: KeyboardEvent) => {
            this.keys.add(e.code);
        });
        document.addEventListener('keyup', (e: KeyboardEvent) => {
            this.keys.delete(e.code);
        });
        this.canvas.addEventListener('click', () => {
            this.canvas.requestPointerLock();
        });
    }

    public tick(): void {
        this.deltaX = 0;
        this.deltaY = 0;
    }
}



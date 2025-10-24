export class InputManager {
    public deltaX = 0;
    public deltaY = 0;
    public keys: Set<string> = new Set();
    public keysPressedThisFrame: Set<string> = new Set();
    public mouseX = 0;
    public mouseY = 0;
    public clicked = false;

    constructor(canvas: HTMLCanvasElement) {
        document.addEventListener('keydown', (e: KeyboardEvent) => {
            // Edge-detect: only mark as pressed-this-frame if it wasn't already down
            if (!this.keys.has(e.code)) {
                this.keysPressedThisFrame.add(e.code);
            }
            this.keys.add(e.code);
        });
        document.addEventListener('keyup', (e: KeyboardEvent) => {
            this.keys.delete(e.code);
        });

        const updateMousePosition = (clientX: number, clientY: number) => {
            const rect = canvas.getBoundingClientRect();
            this.mouseX = clientX - rect.left;
            this.mouseY = clientY - rect.top;
        };

        canvas.addEventListener('pointermove', (e: PointerEvent) => {
            updateMousePosition(e.clientX, e.clientY);
        });
        canvas.addEventListener('pointerdown', (e: PointerEvent) => {
            updateMousePosition(e.clientX, e.clientY);
            this.clicked = true;
        });
    }

    public wasPressed(code: string): boolean {
        return this.keysPressedThisFrame.has(code);
    }

    public tick(): void {
        this.deltaX = 0;
        this.deltaY = 0;
        this.clicked = false;
        this.keysPressedThisFrame.clear();
    }
}



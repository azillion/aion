import "@client/style.css";
import { Renderer } from "@client/renderer/index";
import { UI } from "@client/ui";
import { AppState, CameraMode } from "@client/state";
import { InputManager } from "@client/input";
import { App } from "@client/app";
import { HUDManager } from "@client/hud";
import { CameraManager } from "@client/camera/manager";
import { ClientAuthority } from "@client/authority/clientAuthority";
import { WebWorkerAuthorityProvider, type IAuthorityProvider } from "@client/authority/provider";

// Import all passes and pipelines
import { ShipRelativePipeline } from '@client/renderer/pipelines/shipRelativePipeline';
import { SystemMapPipeline } from '@client/renderer/pipelines/systemMapPipeline';
import type { IRenderPipeline } from '@client/renderer/pipelines/base';

// --- Phase 0: WASM Test ---
async function testWasm() {
    try {
        const response = await fetch('/packages/game-sim/zig-out/bin/game-sim.wasm');
        const wasmBytes = await response.arrayBuffer();
        const { instance } = await WebAssembly.instantiate(wasmBytes);
        const { add } = instance.exports as { add: (a: number, b: number) => number };
        const result = add(2, 3);
        console.log(`[SUCCESS] WASM test passed: 2 + 3 = ${result}`);
        if (result !== 5) {
            throw new Error("WASM add function returned incorrect result!");
        }
    } catch (err) {
        console.error(`[FAILURE] WASM test failed:`, err);
        alert("Failed to load or execute the Zig WASM module. Check the console for errors. You may need to run 'npm run build:sim' first.");
    }
}
// --- End WASM Test ---

async function main() {
    await testWasm(); // Run the test before starting the app
    const canvas = document.getElementById("webgpu-canvas") as HTMLCanvasElement;
    const hudCanvas = document.getElementById("hud-canvas") as HTMLCanvasElement;

    // --- Environment-Agnostic Authority Setup ---
    const authorityProvider: IAuthorityProvider = new WebWorkerAuthorityProvider();
    const connector = await authorityProvider.createConnection();
    const authority = new ClientAuthority(connector);
    // --- END ---

    // --- Application State ---
    const state = new AppState();
    state.playerShipId = 'player-ship';

    // --- Composition Root ---
    
    // 1. Instantiate all rendering pipelines
    const pipelines: Record<CameraMode, IRenderPipeline> = {
      [CameraMode.SHIP_RELATIVE]: new ShipRelativePipeline(),
      [CameraMode.SYSTEM_MAP]: new SystemMapPipeline(),
    };

    // 3. Instantiate core modules, injecting dependencies
    const input = new InputManager(canvas);
    const hud = new HUDManager(hudCanvas);
    const renderer = new Renderer(canvas, state, hud, pipelines);
    const cameraManager = new CameraManager(canvas);
    const ui = new UI(renderer, state, authority, cameraManager);
    
    renderer.ui = ui;

    const app = new App(renderer, authority, state, input, ui, hud, cameraManager);
    await app.start();
}

main().catch(e => {
    console.error(e);
    alert(e);
});

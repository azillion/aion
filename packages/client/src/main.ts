import "@client/style.css";
import { Renderer } from "@client/renderer/index";
import { UI } from "@client/ui";
import { AppState, CameraMode } from "@client/state";
import { InputManager } from "@client/input";
import { App } from "@client/app";
import { HUDManager } from "@client/hud";
import { CameraManager } from "@client/camera/manager";
import { ClientAuthority } from "@client/authority/clientAuthority";
import { WasmAuthorityProvider } from "@client/authority/wasmProvider";
import type { IAuthorityProvider } from "@client/authority/provider";

// Import all passes and pipelines
import { ShipRelativePipeline } from '@client/renderer/pipelines/shipRelativePipeline';
import { SystemMapPipeline } from '@client/renderer/pipelines/systemMapPipeline';
import type { IRenderPipeline } from '@client/renderer/pipelines/base';

declare global {
    interface Window { debugConnector: any }
}

async function main() {
    const canvas = document.getElementById("webgpu-canvas") as HTMLCanvasElement;
    const hudCanvas = document.getElementById("hud-canvas") as HTMLCanvasElement;

    // --- Environment-Agnostic Authority Setup ---
    const authorityProvider: IAuthorityProvider = new WasmAuthorityProvider();
    const connector = await authorityProvider.createConnection();
    const authority = new ClientAuthority(connector);
    // --- END ---

    // --- EXPOSE FOR DEBUGGING ---
    (window as any).debugConnector = connector;
    console.log('Debug connector exposed as `window.debugConnector`. Call `postMessage({ type: \u0064ebugPrintState })` to log Zig state.');
    // --- END DEBUGGING ---

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

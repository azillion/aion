import "./style.css";
import { Renderer } from "./renderer/index";
import { LocalAuthority } from "@server/local";
import { UI } from "./ui";
import { AppState } from "./state";
import { InputManager } from "./input";
import { App } from "./app";
import { HUDManager } from "./hud";
import { CameraManager } from "./camera/manager";

async function main() {
    const canvas = document.getElementById("webgpu-canvas") as HTMLCanvasElement;
    const hudCanvas = document.getElementById("hud-canvas") as HTMLCanvasElement;

    // Instantiate all modules
    const authority = new LocalAuthority();
    const state = new AppState();
    state.playerShipId = 'player-ship';
    const input = new InputManager(canvas);
    const hud = new HUDManager(hudCanvas);
    const renderer = new Renderer(canvas, state, hud);
    const cameraManager = new CameraManager(canvas);
    const ui = new UI(renderer, state, authority, cameraManager);
    
    // The renderer needs a reference to the UI for now
    renderer.ui = ui;

    // Create and start the main App
    const app = new App(renderer, authority, state, input, ui, hud, cameraManager);
    await app.start();
}

main().catch(e => {
    console.error(e);
    alert(e);
});

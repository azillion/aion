import "./style.css";
import { Renderer } from "./renderer";
import { LocalAuthority } from "./authority/local";
import { UI } from "./ui";
import { AppState } from "./state";
import { InputManager } from "./input";

async function main() {
    const canvas = document.getElementById("webgpu-canvas") as HTMLCanvasElement;
    const authority = new LocalAuthority();
    const state = new AppState();
    state.playerShipId = 'player-ship';
    const input = new InputManager(canvas);
    const renderer = new Renderer(canvas, authority, state, input);
    new UI(renderer, state);
    await renderer.start();
}

main().catch(e => {
    console.error(e);
    alert(e);
});

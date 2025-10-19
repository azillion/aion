import "./style.css";
import { Renderer } from "./renderer";
let renderTime = performance.now();

// initWebGPU moved into Renderer

async function main() {
    const canvas = document.getElementById("canvas") as HTMLCanvasElement;
    const renderer = new Renderer(canvas);
    await renderer.start();
}

main().catch(e => {
    console.error(e);
    alert(e);
});

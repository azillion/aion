import "./style.css";
import { Renderer } from "./renderer";
import { LocalAuthority } from "./authority/local";
let renderTime = performance.now();

// initWebGPU moved into Renderer

async function main() {
    const canvas = document.getElementById("canvas") as HTMLCanvasElement;
    const authority = new LocalAuthority();
    const renderer = new Renderer(canvas, authority);
    await renderer.start();
}

main().catch(e => {
    console.error(e);
    alert(e);
});

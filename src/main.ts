import "./style.css";
import { Renderer } from "./renderer";
import { LocalAuthority } from "./authority/local";
import { UI } from "./ui";

// initWebGPU moved into Renderer

async function main() {
    const canvas = document.getElementById("canvas") as HTMLCanvasElement;
    const authority = new LocalAuthority();
    const renderer = new Renderer(canvas, authority);
    new UI(authority);
    await renderer.start();
}

main().catch(e => {
    console.error(e);
    alert(e);
});

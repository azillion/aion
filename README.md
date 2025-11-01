# Project AION (Working Title)

Project AION is an in-development, galaxy-scale, persistent universe simulation built from the ground up with modern web technologies. The project's core is a high-performance simulation engine designed for creating deeply detailed, living worlds with dynamic climates, resource flows, and evolving populations.

The entire simulation is architected for a deep **single-player-first** experience, running locally with maximum performance, while retaining a core design that enables future multiplayer capabilities on an authoritative server.

---

## Core Technical Pillars

This project is an exploration of modern, high-performance game development techniques, built on a unique technology stack:

1.  **Dual-Grid, Hybrid Simulation Engine:** To achieve both deep simulation realism and high-performance rendering, the world is represented on two distinct grids:
    *   **CPU "Simulation Grid":** A coarse-resolution **Goldberg Polyhedron** (or HEALPix-like) grid provides an isotropic environment perfect for running complex, interdependent simulations. This is where deep, lower-frequency systems like global climate, plate tectonics stand-ins, large-scale hydrology, and agent-based economics (faction AI, trade) are modeled with high fidelity.
    *   **GPU "Rendering Grid":** An extremely high-resolution **Cube Sphere** provides the canvas for what the player sees. Its quad-based structure is GPU-native, mapping perfectly to textures and enabling hierarchical Level of Detail (LOD) systems like Quadtree Clipmaps, essential for scaling from orbit to ground level.

2.  **Unified Simulation Core (Zig + WGPU):** The entire simulation logic is written once in **Zig**, a modern, high-performance systems programming language. This single codebase is then:
    *   Compiled to **WebAssembly (WASM)** to run client-side in a Web Worker for seamless, high-performance single-player sessions.
    *   Compiled to a **native library** to serve as the core for future dedicated server applications.

3.  **Real-Time Ray-Traced Rendering (WebGPU):** The client is a sophisticated rendering engine built with the WebGPU API. Instead of traditional rasterization, it utilizes **real-time ray tracing** (e.g., ray marching implicit surfaces/SDFs) to render the world. This approach is a natural fit for procedural generation, allowing for:
    *   **Infinite Detail:** Terrain is represented as a mathematical function, not a triangle mesh, enabling seamless and continuous LOD from orbit to street-level.
    *   **Clean Object Integration:** Procedurally generated structures, cities, and other objects can be integrated perfectly into the world without the complex mesh-stitching or Z-fighting issues of rasterization.

4.  **The Simulation Bridge:** A sophisticated two-way data bridge connects the simulation and rendering grids:
    *   **Physics → Render:** The renderer generates its high-resolution world by performing **barycentric interpolation** on the coarse simulation grid. This smoothly translates the deep simulation data (e.g., temperature, biome type, base elevation) into a continuous field that informs the final render.
    *   **Render → Physics:** Player interactions and other high-frequency events on the fine grid are pushed back to the coarse simulation grid using **area-weighted scattering**, ensuring that player actions have a persistent, meaningful impact on the world state.

5.  **Galaxy-Scale Architecture:** The engine is designed with a tiered simulation lifecycle (`Active`, `Ticking`, `Dormant`) to manage a universe of planets efficiently. This allows the galaxy to feel persistent and alive, with systems evolving even when no players are present, without requiring an impossible amount of computational resources.

---

## Project Structure

The project is a monorepo composed of several key packages:

*   `packages/game-sim/`: **(Zig)** The core of the universe. A self-contained Zig library that holds all simulation logic for both the CPU and GPU layers. It compiles to both WASM and native targets.
*   `packages/client/`: **(TypeScript)** The browser-based client. It handles rendering (WebGPU), user input, and communication with the local simulation host (the WASM module running in a Web Worker).
*   `packages/server/`: **(TypeScript/Bun/Node)** A package reserved for hosting the natively compiled Zig simulation for future multiplayer experiences. In the single-player model, its role is deferred.
*   `packages/shared/`: A shared package defining the data contract (types, constants, network messages) between all other packages.

---

## Current Status: Foundational Prototype Complete

The foundational architecture is now proven and operational. The current prototype successfully demonstrates:
*   End-to-end compilation of the Zig simulation core to both WASM and native.
*   A TypeScript client that loads the WASM module in a Web Worker.
*   The WASM module driving a N-body physics simulation for a solar system, including player ship controls.
*   A multi-scale 3D planet renderer (WebGPU) that visualizes the simulation state in real-time with procedural atmospheric scattering.

## Next Steps

With the core pipeline established, development is focused on realizing the dual-grid architecture and building out the deep planetary simulation.

*   **Phase 4: Planetary Simulation & Rendering:**
    *   **Sim Grid:** Implement the coarse Goldberg/HEALPix grid and the foundational, CPU-based procedural generation step (bedrock, climate, initial soils).
    *   **Render Grid:** Build out the GPU-based ray-marching renderer for the cube-sphere grid.
    *   **The Bridge:** Develop the barycentric interpolation and scattering systems to link the two grids.
    *   **Hydraulics:** Implement the first layer of the real-time GPU substrate simulation, focusing on hydraulic flow and erosion.

*   **Phase 5:** Implementing the CPU "Agent Layer" and enabling two-way interaction, allowing player actions to have a persistent effect on the world.

*   **Phase 6:** Building out the server-side lifecycle management (`Active`/`Ticking`/`Dormant` states) to enable galaxy-scale persistence in the single-player universe.
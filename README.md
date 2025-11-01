# Project AION (Working Title)

Project AION is an in-development, galaxy-scale, persistent universe simulation built from the ground up with modern web technologies. The project's core is a high-performance, GPU-accelerated simulation engine designed for creating detailed, living worlds with dynamic climates, resource flows, and evolving populations.

The entire simulation is architected to be multiplayer-first, running on an authoritative server, with a seamless transition for local single-player experiences.

---

## Core Technical Pillars

This project is an exploration of modern, high-performance game development techniques, built on a unique technology stack:

1.  **Data-Oriented, Hybrid Simulation Engine:** The simulation is split into two layers for maximum performance and complexity:
    *   **GPU "Substrate" Simulation:** A massively parallel cellular automaton runs directly on the GPU, simulating high-frequency, local physical phenomena like water flow, erosion, and heat diffusion across millions of cells on a planet's surface.
    *   **CPU "Agent" Simulation:** A low-frequency, high-level simulation written in Zig manages global, intelligent systems like economics, faction AI, and trade routes, interacting with the GPU substrate through a clean command interface.

2.  **Unified Simulation Core (Zig + WGPU):** The entire simulation logic is written once in **Zig**, a modern, high-performance systems programming language. This single codebase is then:
    *   Compiled to **WebAssembly (WASM)** to run client-side in a Web Worker for seamless, high-performance local/single-player sessions.
    *   Compiled to a **native library** to run at maximum speed on a dedicated server for the authoritative multiplayer experience.

3.  **Procedural & Stateful Hybrid Model:** To manage memory efficiently, the engine only stores what cannot be deterministically recalculated. Base terrain and climate are generated procedurally from a seed, while dynamic state (like water levels, population, and player actions) is stored in GPU textures. This allows for incredibly detailed worlds with a minimal memory footprint.

4.  **Galaxy-Scale Architecture:** The server is designed with a tiered simulation lifecycle (`Active`, `Ticking`, `Dormant`) to manage a universe of planets. This allows the galaxy to feel persistent and alive, with systems evolving even when no players are present, without requiring an impossible amount of server resources.

5.  **Modern Web Rendering (WebGPU):** The client is a sophisticated rendering engine built with the WebGPU API, designed to visualize the complex state of the simulation in real-time with procedural terrain, atmospheric scattering, and dynamic surface features driven by the simulation state.

---

## Project Structure

The project is a monorepo composed of several key packages:

*   `packages/game-sim/`: **(Zig)** The core of the universe. A self-contained Zig library that holds all simulation logic for both the GPU and CPU layers. It compiles to both WASM and native targets.
*   `packages/client/`: **(TypeScript)** The browser-based client. It handles rendering (WebGPU), user input, and communication with a simulation host. For local play, the host is the WASM module running in a Web Worker.
*   `packages/server/`: **(TypeScript/Bun/Node)** The authoritative server. It manages the world state, player connections, and hosts the natively compiled Zig simulation library for multiplayer games.
*   `packages/shared/`: A shared package defining the data contract (types, constants, network messages) between all other packages.

---

## Current Status: Foundational Prototype Complete

The foundational architecture is now proven and operational. The current prototype successfully demonstrates:
*   End-to-end compilation of the Zig simulation core to both WASM and native.
*   A TypeScript client that loads the WASM module in a Web Worker.
*   The WASM module driving a N-body physics simulation for a solar system, including player ship controls.
*   A multi-scale 3D planet renderer (WebGPU) that visualizes the simulation state in real-time with procedural atmospheric scattering.

## Next Steps

With the core pipeline established, development is focused on expanding the simulation's complexity and adding interactive systems.
*   **Phase 4:** Expand the physical simulation with procedural terrain and hydraulic flow.
*   **Phase 5:** Implementing the CPU "Agent Layer" and the two-way bridge, allowing player interaction to have a persistent effect on the world.
*   **Phase 6:** Building out the server-side lifecycle management (`Active`/`Dormant` states) to enable galaxy-scale persistence.

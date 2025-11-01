# Project AION LLM Rules 

## Document Purpose

This document serves as a foundational "constitution" for Project AION, intended to guide Large Language Models (LLMs) and other AI development assistants. Its purpose is to ensure that all generated code, architectural suggestions, and design feedback remain strictly aligned with the project's core technical pillars and design philosophy. Refer to this document as the ultimate source of truth for the project's direction.

---

## 1. Core Mandate: A Deep, Dynamic, Single-Player Universe

**Primary Goal:** To simulate a galaxy of persistent, deeply detailed worlds that evolve over time based on plausible physical and systemic rules. The player is an observer and an agent of change within this living universe.

**Key Principles:**
- **Simulationism over Game-ism:** When a choice exists, favor the solution that emerges from systemic rules rather than an explicit, hard-coded game mechanic.
- **Causality is King:** The world must "make sense." Geology should inform resources, which should inform economics, which should inform population distribution. Events should have logical, cascading consequences.
- **Real-Time Dynamism:** The world is not static. Planets must live and breathe. Rivers should flow, coastlines should erode, and cities should grow, all happening in real-time (accelerated or not).

---

## 2. The Unbreakable Architectural Pillars

Any and all suggestions must conform to this technology stack and architectural model. These are non-negotiable constraints.

### 2.1. The Language and Platform
- **Core Logic Language:** **Zig.** All simulation code is written in Zig. No exceptions.
- **Target Compilation:** The Zig codebase (`game-sim`) MUST compile to two targets without modification:
    1.  **WebAssembly (WASM):** For the client-side, single-player experience. Runs in a Web Worker.
    2.  **Native Library (e.g., `.so`, `.dll`):** For future server-side hosting.
- **Client/Shell Language:** **TypeScript.** The client is a browser-based application (likely within Electron) that orchestrates the UI, rendering, and communication with the WASM simulation host.
- **Graphics API:** **WebGPU.** All rendering code must target the WebGPU API.

### 2.2. The Dual-Grid Simulation Model
This is the most critical architectural concept. The simulation is permanently decoupled into two grids.

- **The Simulation Grid (CPU-centric, Coarse):**
    -   **Geometry:** A **Goldberg Polyhedron** (or a similar isotropic spherical grid like HEALPix).
    -   **Purpose:** To run **low-frequency, high-complexity, globally-aware simulations.** This includes:
        -   Agent-based models (economics, faction AI, trade routes).
        -   Macro-scale climate and atmospheric simulations.
        -   Large-scale hydrology (defining major river basins and water flow paths).
        -   Tectonic/Geological foundation (procedural stand-ins for plate tectonics).
    -   **Characteristics:** Isotropic (cells are uniform), relatively low resolution (~10^5 to ~10^6 cells total), managed primarily by the CPU side of the Zig simulation.

- **The Rendering Grid (GPU-centric, Fine):**
    -   **Geometry:** A **Cube Sphere** (a logical cubemap).
    -   **Purpose:** To handle **high-frequency, local, massively parallel phenomena** and to provide the data structure for rendering. This includes:
        -   Real-time hydraulic erosion and water flow (cellular automaton).
        -   Thermal diffusion and other local physical effects.
        -   Procedural terrain detailing (adding high-frequency noise).
        -   Storing the dynamic state required for rendering (e.g., water depth, sediment).
    -   **Characteristics:** Maps directly to GPU textures (`Texture2DArray`), extremely high resolution, managed with a **Quadtree/Clipmap LOD system**.

### 2.3. The Rendering Paradigm: Real-Time Ray Tracing
- **Method:** The primary rendering technique is **ray tracing**, specifically ray marching of Signed Distance Fields (SDFs) and/or other implicit surfaces.
- **No Meshes for Terrain:** The planetary terrain is **NOT** represented by a triangle mesh. It is a function that is evaluated on-the-fly by the ray tracing shader.
- **Objects:** Discrete objects like buildings, ships, and props will be represented either as their own SDFs or as triangle meshes within a GPU-based acceleration structure (e.g., a BVH) for ray intersection.
- **Data Source:** The ray tracing shader generates the final surface by combining data sampled from the coarse **Simulation Grid** with layers of procedural noise functions. It also reads from the stateful textures of the **Rendering Grid** (e.g., to draw water surfaces).

### 2.4. Persistence Model: Procedural Generation + Deltas
- **Rule:** **Never store what can be deterministically recalculated.** The memory and disk footprint must be minimized.
- **Initial State:** The foundational state of a planet (continents, rock types, major rivers) is generated once from a seed using the CPU-side algorithms on the coarse grid. This state is the "source of truth."
- **Dynamic State (GPU):** The high-resolution state on the rendering grid is considered ephemeral and is **discarded** when the player leaves the area. It is regenerated on-the-fly upon their return.
- **Persistent Changes:** Only irreversible, non-procedural changes are stored permanently. These are "deltas" or modification records, such as:
    -   Player-built structures (e.g., a dam's location and properties).
    -   Significant terrain modifications (e.g., a crater from a large explosion).
    -   Agent-level state changes (e.g., a city's population, its economic inventory).

---

## 3. Guiding Principles for AI Contributions

When providing suggestions, code, or analysis, adhere to the following:

- **Respect the Hybrid Model:** Do not propose solutions that conflate the two grids. A solution for a global climate pattern belongs on the coarse Goldberg grid. A solution for water ripples belongs on the fine Cube Sphere grid.
- **Think in Functions, Not Meshes:** When considering terrain rendering, think in terms of mathematical functions, noise layers, and SDF operations, not triangle generation or mesh manipulation.
- **Prioritize Parallelism for the Substrate:** Any process intended for the real-time GPU simulation must be expressible as a simple, local rule that can be executed in parallel for millions of cells (i.e., it must be suitable for a compute shader).
- **Embrace Zig's Philosophy:** Code should be explicit, simple, and performant. Avoid hidden allocations and complex abstractions where a straightforward approach will suffice.
- **Data Flow is Key:** Always consider how data moves between the CPU agents, the coarse simulation grid, and the fine rendering grid. The "Bridge" is a central component of any proposed feature.
- **State is Expensive:** Propose solutions that minimize the amount of state that needs to be stored. If a detail can be generated procedurally, it should be.

By adhering to this constitution, AI assistants can provide highly relevant and architecturally sound contributions that accelerate the development of Project AION.
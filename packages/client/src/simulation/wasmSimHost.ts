import type { ISimulationHost, SimulationHandle } from "@shared/simulation";

export class WasmSimHost implements ISimulationHost {
  private simState: SimulationHandle | null = null;
  private _wasmInstance: WebAssembly.Instance | null = null;
  private _create_simulator:
    | ((device: GPUDevice, queue: GPUQueue) => SimulationHandle)
    | null = null;
  private _tick_simulator: ((state: SimulationHandle) => void) | null = null;
  private _get_output_texture_view: ((state: SimulationHandle) => any) | null = null; // We'll handle the return type later
  private _destroy_simulator: ((state: SimulationHandle) => void) | null = null;

  async initialize(device: GPUDevice, queue: GPUQueue): Promise<void> {
    const response = await fetch(
      "/packages/game-sim/zig-out/lib/game-sim.wasm"
    );
    if (!response.ok) {
      throw new Error(`Failed to fetch WASM: ${response.status} ${response.statusText}`);
    }

    const wasmBytes = await response.arrayBuffer();
    const { instance } = await WebAssembly.instantiate(wasmBytes, {});
    this._wasmInstance = instance;

    const exports = instance.exports as Record<string, unknown>;
    this._create_simulator = exports["create_simulator"] as unknown as (
      device: GPUDevice,
      queue: GPUQueue
    ) => SimulationHandle;
    this._tick_simulator = exports["tick_simulator"] as unknown as (
      state: SimulationHandle
    ) => void;
    this._get_output_texture_view = exports["get_output_texture_view"] as unknown as (
      state: SimulationHandle
    ) => any;
    this._destroy_simulator = exports["destroy_simulator"] as unknown as (
      state: SimulationHandle
    ) => void;

    if (!this._create_simulator) {
      throw new Error("WASM export 'create_simulator' not found");
    }

    const simState = this._create_simulator(device, queue);
    this.simState = simState;
    if (!this.simState || this.simState === 0) {
      throw new Error("Simulator failed to initialize (null/zero handle)");
    }

    console.log("WASM simulator initialized", { handle: this.simState });
  }

  tick(): void {
    if (!this.simState) {
      console.warn("tick() called before simulator is initialized");
      return;
    }
    if (!this._tick_simulator) {
      console.warn("tick() requested but WASM export 'tick_simulator' is missing");
      return;
    }
    this._tick_simulator(this.simState);
  }

  getOutputTextureView(): any {
    // TODO: This currently returns a numeric handle from Zig. We'll bridge this
    // to an actual GPUTextureView on the TS side in a subsequent step.
    if (!this.simState) {
      console.warn("getOutputTextureView() called before simulator is initialized");
      return null;
    }
    if (!this._get_output_texture_view) {
      console.warn(
        "getOutputTextureView() requested but WASM export 'get_output_texture_view' is missing"
      );
      return null;
    }
    const handle = this._get_output_texture_view(this.simState);
    return handle;
  }
}



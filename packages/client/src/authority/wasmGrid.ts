type WasmGridExports = {
  memory: WebAssembly.Memory;
  get_query_scratch_ptr: () => number;
  get_scratch_buffer_ptr: () => number;
  create_grid: (size: number) => void;
  get_grid_vertex_buffer_out: (outPtr: number) => void;
  get_grid_elevation_buffer_out: (outPtr: number) => void;
  get_grid_index_buffer_out: (outPtr: number) => void;
};

export class WasmGridBridge {
  private instance: WebAssembly.Instance | null = null;
  private exports!: WasmGridExports;
  private scratchPtr: number = 0;

  public static async create(): Promise<WasmGridBridge> {
    const bridge = new WasmGridBridge();
    await bridge.initialize();
    return bridge;
  }

  private async initialize(): Promise<void> {
    const response = await fetch('/packages/game-sim/zig-out/bin/game-sim.wasm');
    const bytes = await response.arrayBuffer();
    const importObject = {
      env: {
        log_error: (ptr: number, len: number) => {
          if (!this.instance) return;
          const mem = (this.instance.exports as any).memory as WebAssembly.Memory;
          const view = new Uint8Array(mem.buffer, ptr, len);
          const msg = new TextDecoder().decode(view);
          console.error(`[Zig WASM Error]: ${msg}`);
        },
      }
    } as any;
    const { instance } = await WebAssembly.instantiate(bytes, importObject);
    this.instance = instance;
    this.exports = instance.exports as unknown as WasmGridExports;
    this.scratchPtr = this.exports.get_query_scratch_ptr();
  }

  public createGrid(size: number): void {
    this.exports.create_grid(size >>> 0);
  }

  private readSlice(outPtr: number): { ptr: number; len: number } {
    const dv = new DataView(this.exports.memory.buffer);
    const ptr = dv.getUint32(outPtr + 0, true);
    const len = dv.getUint32(outPtr + 4, true);
    return { ptr, len };
  }

  public getGridVertexBuffer(): Float32Array {
    this.exports.get_grid_vertex_buffer_out(this.scratchPtr);
    const { ptr, len } = this.readSlice(this.scratchPtr);
    if (!ptr || !len) return new Float32Array(0);
    return new Float32Array(this.exports.memory.buffer, ptr, len / 4);
  }

  public getGridElevationBuffer(): Float32Array {
    this.exports.get_grid_elevation_buffer_out(this.scratchPtr);
    const { ptr, len } = this.readSlice(this.scratchPtr);
    if (!ptr || !len) return new Float32Array(0);
    return new Float32Array(this.exports.memory.buffer, ptr, len / 4);
  }

  public getGridIndexBuffer(): Uint32Array {
    this.exports.get_grid_index_buffer_out(this.scratchPtr);
    const { ptr, len } = this.readSlice(this.scratchPtr);
    if (!ptr || !len) return new Uint32Array(0);
    return new Uint32Array(this.exports.memory.buffer, ptr, len / 4);
  }
}



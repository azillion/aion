type WasmGridExports = {
  memory: WebAssembly.Memory;
  get_scratch_buffer_ptr: () => number;
  get_scratch_buffer_len: () => number;
  load_prebaked_planet: (ptr: number, len: number) => boolean;
  unload_prebaked_planet: () => void;
  get_planet_resolution: () => number;
  get_planet_tile_count: () => number;
  get_planet_tile_stride: () => number;
  read_planet_tile: (tileIndex: number, outPtr: number) => boolean;
  get_planet_neighbor: (tileIndex: number, dir: number) => number;
  get_planet_vertex_count: () => number;
  get_planet_vertices_ptr: () => number;
  get_planet_elevations_ptr: () => number;
  get_planet_index_count: () => number;
  get_planet_indices_ptr: () => number;
  get_planet_edge_index_count: () => number;
  get_planet_edge_indices_ptr: () => number;
};

export type PlanetBuffers = {
  vertices: Float32Array;
  elevations: Float32Array;
  indices: Uint32Array;
  edges: Uint32Array;
  resolution: number;
};

export class WasmGridBridge {
  private instance: WebAssembly.Instance | null = null;
  private exports!: WasmGridExports;
  private scratchPtr = 0;
  private scratchLen = 0;

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
    this.scratchPtr = this.exports.get_scratch_buffer_ptr();
    this.scratchLen = this.exports.get_scratch_buffer_len();
  }

  private getMemoryU8(): Uint8Array {
    return new Uint8Array(this.exports.memory.buffer);
  }

  public async loadPlanetFromUrl(url: string): Promise<PlanetBuffers> {
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`Failed to fetch planet blob: ${response.status} ${response.statusText}`);
    }
    const buffer = await response.arrayBuffer();
    const bytes = new Uint8Array(buffer);
    if (bytes.length > this.scratchLen) {
      throw new Error(`Planet blob (${bytes.length} bytes) exceeds scratch buffer (${this.scratchLen} bytes)`);
    }

    this.getMemoryU8().set(bytes, this.scratchPtr);
    if (!this.exports.load_prebaked_planet(this.scratchPtr, bytes.length)) {
      throw new Error('load_prebaked_planet returned false');
    }

    const vertexCount = this.exports.get_planet_vertex_count();
    const vertexPtr = this.exports.get_planet_vertices_ptr();
    const elevationPtr = this.exports.get_planet_elevations_ptr();
    const indexCount = this.exports.get_planet_index_count();
    const indexPtr = this.exports.get_planet_indices_ptr();
    const edgeCount = this.exports.get_planet_edge_index_count();
    const edgePtr = this.exports.get_planet_edge_indices_ptr();

    const vertices = vertexCount > 0 && vertexPtr !== 0
      ? new Float32Array(new Float32Array(this.exports.memory.buffer, vertexPtr, vertexCount * 3))
      : new Float32Array();
    const elevations = vertexCount > 0 && elevationPtr !== 0
      ? new Float32Array(new Float32Array(this.exports.memory.buffer, elevationPtr, vertexCount))
      : new Float32Array();
    const indices = indexCount > 0 && indexPtr !== 0
      ? new Uint32Array(new Uint32Array(this.exports.memory.buffer, indexPtr, indexCount))
      : new Uint32Array();
    const edges = edgeCount > 0 && edgePtr !== 0
      ? new Uint32Array(new Uint32Array(this.exports.memory.buffer, edgePtr, edgeCount))
      : new Uint32Array();

    return {
      vertices,
      elevations,
      indices,
      edges,
      resolution: this.exports.get_planet_resolution(),
    };
  }

  public dispose(): void {
    this.exports.unload_prebaked_planet();
  }
}



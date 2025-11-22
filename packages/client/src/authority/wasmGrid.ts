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
};

export type PlanetBuffers = {
  vertices: Float32Array;
  elevations: Float32Array;
  indices: Uint32Array;
  resolution: number;
};

const UNSET_NEIGHBOR = 0xFFFFFFFF;
const DISK_TILE_OFFSETS = {
  pos: 0,
  face: 12,
  axialQ: 16,
  axialR: 20,
  neighbors: 24,
  flags: 48,
} as const;

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

  private getScratchView(byteLength: number): DataView {
    return new DataView(this.exports.memory.buffer, this.scratchPtr, byteLength);
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

    const tileCount = this.exports.get_planet_tile_count();
    const stride = this.exports.get_planet_tile_stride();
    const vertices = new Float32Array(tileCount * 3);
    const elevations = new Float32Array(tileCount); // Placeholder until we bake height data offline
    const triangleList: number[] = [];
    const neighborScratch = new Array<number>(6);

    for (let i = 0; i < tileCount; i++) {
      if (!this.exports.read_planet_tile(i, this.scratchPtr)) {
        throw new Error(`read_planet_tile failed for tile ${i}`);
      }
      const view = this.getScratchView(stride);
      const base = i * 3;
      vertices[base + 0] = view.getFloat32(DISK_TILE_OFFSETS.pos + 0, true);
      vertices[base + 1] = view.getFloat32(DISK_TILE_OFFSETS.pos + 4, true);
      vertices[base + 2] = view.getFloat32(DISK_TILE_OFFSETS.pos + 8, true);

      let neighborCount = 0;
      for (let dir = 0; dir < 6; dir++) {
        const nb = view.getUint32(DISK_TILE_OFFSETS.neighbors + dir * 4, true);
        if (nb !== UNSET_NEIGHBOR) {
          neighborScratch[neighborCount++] = nb;
        }
      }
      this.appendTrianglesForTile(i, neighborScratch, neighborCount, triangleList);
    }

    return {
      vertices,
      elevations,
      indices: new Uint32Array(triangleList),
      resolution: this.exports.get_planet_resolution(),
    };
  }

  private appendTrianglesForTile(tileIndex: number, neighbors: number[], neighborCount: number, out: number[]): void {
    if (neighborCount < 2) return;
    for (let i = 0; i < neighborCount; i++) {
      const b = neighbors[i];
      const c = neighbors[(i + 1) % neighborCount];
      if (tileIndex < b && tileIndex < c) {
        out.push(tileIndex, b, c);
      }
    }
  }

  public dispose(): void {
    this.exports.unload_prebaked_planet();
  }
}



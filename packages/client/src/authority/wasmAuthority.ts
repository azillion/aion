import type { Authority, InputState } from '@shared/authority';
import type { SystemState } from '@shared/types';

type WasmExports = {
    memory: WebAssembly.Memory;
    create_simulator: (ptr: number, len: number) => number;
    tick_simulator: (simPtr: number, dt: number, inputPtr: number) => void;
    query_simulator_state: (simPtr: number) => bigint;
    free_result_buffer: (addr: number, len: number) => void;
    destroy_simulator: (simPtr: number) => void;
};

export class WasmAuthority implements Authority {
    private instance: WebAssembly.Instance | null = null;
    private exports!: WasmExports;
    private simPtr: number = 0;

    constructor() {}

    async initialize(initialState: SystemState): Promise<void> {
        const response = await fetch('/packages/game-sim/zig-out/lib/game-sim.wasm');
        const bytes = await response.arrayBuffer();
        const { instance } = await WebAssembly.instantiate(bytes, {});
        this.instance = instance;
        const e = instance.exports as unknown as WasmExports;
        this.exports = e;

        const json = JSON.stringify(initialState);
        const enc = new TextEncoder();
        const encoded = enc.encode(json);
        const mem = new Uint8Array(e.memory.buffer);
        const base = 100000; // temporary fixed offset; replace with proper allocator later
        mem.set(encoded, base);
        this.simPtr = e.create_simulator(base, encoded.length);
        if (!this.simPtr) throw new Error('Failed to create simulator');
    }

    async query(): Promise<SystemState> {
        if (!this.instance) throw new Error('WASM not initialized');
        const packed = this.exports.query_simulator_state(this.simPtr);
        const addr = Number((packed >> 32n) & 0xffffffffn);
        const len = Number(packed & 0xffffffffn);
        const mem = new Uint8Array(this.exports.memory.buffer, addr, len);
        const json = new TextDecoder().decode(mem);
        this.exports.free_result_buffer(addr, len);
        return JSON.parse(json) as SystemState;
    }

    async tick(deltaTime: number, input: InputState): Promise<void> {
        if (!this.instance) throw new Error('WASM not initialized');
        // Build keys mask (must match Zig bit layout)
        let mask = 0 >>> 0;
        const set = input.keys;
        const bit = (b: number) => (mask |= (1 << b) >>> 0);
        if (set.has('KeyW')) bit(0);
        if (set.has('KeyS')) bit(1);
        if (set.has('KeyA')) bit(2);
        if (set.has('KeyD')) bit(3);
        if (set.has('Space')) bit(4);
        if (set.has('ControlLeft')) bit(5);
        if (set.has('Backspace')) bit(6);
        if (set.has('KeyR')) bit(7);
        if (set.has('KeyF')) bit(8);
        if (set.has('KeyZ')) bit(9);
        if (set.has('KeyC')) bit(10);
        if (set.has('KeyQ')) bit(11);
        if (set.has('KeyE')) bit(12);
        if (set.has('KeyX')) bit(13);
        if (set.has('Backquote')) bit(14);

        const base = 120000; // temporary input struct location
        const mem = new DataView(this.exports.memory.buffer);
        mem.setUint32(base, mask >>> 0, true);
        this.exports.tick_simulator(this.simPtr, deltaTime, base);
    }

    async setTimeScale(_scale: number): Promise<void> {}
    async addBody(_body: Omit<SystemState['bodies'][number], 'id'>): Promise<void> {}
    async autoLand(_targetBodyId: string | null): Promise<void> {}
    async teleportToSurface(_targetBodyId: string | null): Promise<void> {}
}



import type { Authority, InputState } from '@shared/authority';
import type { SystemState } from '@shared/types';

type WasmExports = {
    memory: WebAssembly.Memory;
    get_input_buffer: () => number;
    get_scratch_buffer_ptr: () => number;
    create_simulator: (ptr: number, len: number) => number;
    tick_simulator: (simPtr: number, dt: number) => void;
    // Out-param variant to avoid ABI ambiguity
    query_simulator_state_out: (simPtr: number, outPtr: number) => void;
    free_result_buffer: (addr: number, len: number) => void;
    destroy_simulator: (simPtr: number) => void;
    // Scratch pointer provider from Zig (8 bytes for [ptr,len])
    get_query_scratch_ptr: () => number;
    // Control and actions
    set_time_scale: (simPtr: number, scale: number) => void;
    add_body: (simPtr: number, jsonPtr: number, jsonLen: number) => void;
    teleport_to_surface: (simPtr: number, idPtr: number, idLen: number) => void;
    auto_land: (simPtr: number, idPtr: number, idLen: number) => void;
};

export class WasmAuthority implements Authority {
    private instance: WebAssembly.Instance | null = null;
    private exports!: WasmExports;
    private simPtr: number = 0;
    private scratchPtr: number = 0;
    private scratchBufPtr: number = 0;

    constructor() {}

    async initialize(initialState: SystemState): Promise<void> {
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
            },
        } as any;
        const { instance } = await WebAssembly.instantiate(bytes, importObject);
        this.instance = instance;
        const e = instance.exports as unknown as WasmExports;
        this.exports = e;
        this.scratchPtr = e.get_query_scratch_ptr();
        this.scratchBufPtr = e.get_scratch_buffer_ptr();

        const json = JSON.stringify(initialState);
        console.log(json);
        const encoded = new TextEncoder().encode(json);
        // Write at a fixed offset in linear memory for initialization
        const base = 100000; // temporary fixed offset for init JSON
        new Uint8Array(e.memory.buffer).set(encoded, base);

        this.simPtr = e.create_simulator(base, encoded.length);
        if (!this.simPtr) throw new Error('Failed to create simulator');
    }

    async query(): Promise<SystemState> {
        if (!this.instance) throw new Error('WASM not initialized');
        this.exports.query_simulator_state_out(this.simPtr, this.scratchPtr);
        const dv = new DataView(this.exports.memory.buffer);
        const ptr = dv.getUint32(this.scratchPtr + 0, true);
        const len = dv.getUint32(this.scratchPtr + 4, true);
        if (!ptr || !len) throw new Error('query_simulator_state returned empty buffer');
        const mem = new Uint8Array(this.exports.memory.buffer, ptr, len);
        const copy = new Uint8Array(len);
        copy.set(mem);
        this.exports.free_result_buffer(ptr, len);
        const json = new TextDecoder().decode(copy);
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
        if (set.has('KeyB')) bit(15);

        const inputPtr = this.exports.get_input_buffer();
        new DataView(this.exports.memory.buffer).setUint32(inputPtr, mask >>> 0, true);
        this.exports.tick_simulator(this.simPtr, deltaTime);
    }

    async setTimeScale(scale: number): Promise<void> {
        if (!this.instance) throw new Error('WASM not initialized');
        this.exports.set_time_scale(this.simPtr, scale);
    }

    async addBody(body: Omit<SystemState['bodies'][number], 'id'>): Promise<void> {
        if (!this.instance) throw new Error('WASM not initialized');
        const json = JSON.stringify(body);
        const encoded = new TextEncoder().encode(json);
        new Uint8Array(this.exports.memory.buffer, this.scratchBufPtr, encoded.length).set(encoded);
        this.exports.add_body(this.simPtr, this.scratchBufPtr, encoded.length);
    }

    async autoLand(targetBodyId: string | null): Promise<void> {
        if (!this.instance || !targetBodyId) return;
        const encoded = new TextEncoder().encode(targetBodyId);
        new Uint8Array(this.exports.memory.buffer, this.scratchBufPtr, encoded.length).set(encoded);
        new Uint8Array(this.exports.memory.buffer, this.scratchBufPtr + encoded.length, Math.max(0, 4096 - encoded.length)).fill(0);
        this.exports.auto_land(this.simPtr, this.scratchBufPtr, encoded.length);
    }

    async teleportToSurface(targetBodyId: string | null): Promise<void> {
        if (!this.instance || !targetBodyId) return;
        const encoded = new TextEncoder().encode(targetBodyId);
        new Uint8Array(this.exports.memory.buffer, this.scratchBufPtr, encoded.length).set(encoded);
        new Uint8Array(this.exports.memory.buffer, this.scratchBufPtr + encoded.length, Math.max(0, 4096 - encoded.length)).fill(0);
        this.exports.teleport_to_surface(this.simPtr, this.scratchBufPtr, encoded.length);
    }
}



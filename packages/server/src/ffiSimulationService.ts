import { dlopen, FFIType, suffix, toArrayBuffer } from 'bun:ffi';
import type { Pointer } from 'bun:ffi';
import type { ISimulationService } from './simulationService';

const GRID_WIDTH = 2048;
const GRID_HEIGHT = 2048;
const GRID_SIZE_BYTES = GRID_WIDTH * GRID_HEIGHT * 4;

const libPath = new URL(`../zig-out/lib/libsimulation-service.${suffix}`, import.meta.url).pathname;

export class FfiSimulationService implements ISimulationService {
    private simStatePtr: Pointer | null = null;

    private symbols: {
        create_simulator: () => Pointer;
        tick_simulator: (ptr: Pointer) => void;
        get_simulation_data: (ptr: Pointer) => Pointer;
        release_simulation_data: (ptr: Pointer) => void;
        destroy_simulator: (ptr: Pointer) => void;
    };

    constructor() {
        const { symbols } = dlopen(libPath, {
            create_simulator: { returns: FFIType.pointer },
            tick_simulator: { args: [FFIType.pointer] },
            // treat data pointer as a pointer, not u64
            get_simulation_data: { args: [FFIType.pointer], returns: FFIType.pointer },
            release_simulation_data: { args: [FFIType.pointer] },
            destroy_simulator: { args: [FFIType.pointer] },
        });
        this.symbols = symbols as any;
    }

    async initialize(): Promise<void> {
        this.simStatePtr = this.symbols.create_simulator();
        if (!this.simStatePtr || this.simStatePtr === 0) {
            throw new Error('Failed to initialize Zig simulation service (returned null pointer).');
        }
    }

    tick(): void {
        if (!this.simStatePtr) return;
        this.symbols.tick_simulator(this.simStatePtr);
    }

    getLatestData(): Buffer | null {
        if (!this.simStatePtr) return null;

        const dataPtr = this.symbols.get_simulation_data(this.simStatePtr);
        // const ptrBI = BigInt(dataPtr as unknown as bigint);

        if (!dataPtr || Number(dataPtr) === 0) {
            // eslint-disable-next-line no-console
            console.error('Bun: Data pointer from Zig is null.');
            this.symbols.release_simulation_data(this.simStatePtr);
            return null;
        }
        try {
            // eslint-disable-next-line no-console
            const ab = toArrayBuffer(dataPtr, 0, GRID_SIZE_BYTES);
            const buf = Buffer.from(ab);
            this.symbols.release_simulation_data(this.simStatePtr);

            if (buf.length !== GRID_SIZE_BYTES) {
                // eslint-disable-next-line no-console
                console.error(`wrong size: expected ${GRID_SIZE_BYTES}, got ${buf.length}`);
                return null;
            }
            return buf;
        } catch (e) {
            // eslint-disable-next-line no-console
            console.error('Bun: CRITICAL ERROR creating Uint8Array from FFI pointer.', e);
            this.symbols.release_simulation_data(this.simStatePtr);
            return null;
        }
    }

    shutdown(): void {
        if (this.simStatePtr) {
            this.symbols.destroy_simulator(this.simStatePtr);
            this.simStatePtr = null;
        }
    }
}



const std = @import("std");

const SimulatorState = struct {
    _dummy: u8,
};

fn wasmAllocator() std.mem.Allocator {
    const builtin = @import("builtin");
    return if (builtin.target.cpu.arch == .wasm32) std.heap.wasm_allocator else std.heap.page_allocator;
}

export fn create_simulator(_: usize, _: usize) callconv(.c) ?*SimulatorState {
    const allocator = wasmAllocator();
    const state = allocator.create(SimulatorState) catch return null;
    state.* = .{ ._dummy = 0 };
    return state;
}

export fn tick_simulator(_: *SimulatorState) callconv(.c) void {
    // no-op in CPU-only stub
}

export fn get_output_texture_view(_: *SimulatorState) callconv(.c) usize {
    return 0;
}

export fn destroy_simulator(state: *SimulatorState) callconv(.c) void {
    const allocator = wasmAllocator();
    allocator.destroy(state);
}

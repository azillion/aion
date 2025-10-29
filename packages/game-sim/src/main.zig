const std = @import("std");
const build_options = @import("build_options");
const wgpu = blk: {
    if (build_options.enable_wgpu_headers) {
        break :blk @cImport(@cInclude("webgpu/wgpu.h"));
    } else {
        break :blk struct {};
    }
};

// This `export` keyword makes the function visible to the linker.
export fn add(a: i32, b: i32) i32 {
    return a + b;
}

test "basic add functionality" {
    try std.testing.expect(add(2, 2) == 4);
}

test "wgpu header types available" {
    if (build_options.enable_wgpu_headers) {
        const SURFACE: wgpu.WGPUSurface = undefined;
        _ = SURFACE;
    }
    try std.testing.ok(true);
}

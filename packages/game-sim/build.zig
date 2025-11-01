const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});

    const wasm_target = b.standardTargetOptions(.{ .default_target = .{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    } });

    const wasm_exe = b.addExecutable(.{
        .name = "game-sim",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = wasm_target,
            .optimize = optimize,
        }),
    });

    // Add our new directories as modules (Zig 0.15.x API)
    wasm_exe.root_module.addImport("core", b.createModule(.{ .root_source_file = b.path("src/core/types.zig") }));
    wasm_exe.root_module.addImport("ffi", b.createModule(.{ .root_source_file = b.path("src/ffi/index.zig") }));
    wasm_exe.root_module.addImport("math", b.createModule(.{ .root_source_file = b.path("src/math/index.zig") }));
    wasm_exe.root_module.addImport("physics", b.createModule(.{ .root_source_file = b.path("src/physics/nbody.zig") }));
    wasm_exe.root_module.addImport("sim", b.createModule(.{ .root_source_file = b.path("src/sim/index.zig") }));

    // Ensure required FFI symbols are exported and not stripped
    wasm_exe.rdynamic = true;

    b.installArtifact(wasm_exe);
}

const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    });

    // root module for the exe is the ffi entry
    const root_mod = b.createModule(.{
        .root_source_file = b.path("src/ffi/index.zig"),
        .target = target,
        .optimize = optimize,
    });

    const wasm_exe = b.addExecutable(.{
        .name = "game-sim",
        .root_module = root_mod,
    });

    // deps as modules (with target/optimize)
    const math_mod = b.createModule(.{ .root_source_file = b.path("src/math/index.zig"), .target = target, .optimize = optimize });
    const core_mod = b.createModule(.{ .root_source_file = b.path("src/core/types.zig"), .target = target, .optimize = optimize });
    const physics_mod = b.createModule(.{ .root_source_file = b.path("src/physics/nbody.zig"), .target = target, .optimize = optimize });
    const sim_mod = b.createModule(.{ .root_source_file = b.path("src/sim/index.zig"), .target = target, .optimize = optimize });
    const planet_mod = b.createModule(.{ .root_source_file = b.path("src/planet/index.zig"), .target = target, .optimize = optimize });

    // wire module dependencies transitively
    // core -> math
    core_mod.addImport("math", math_mod);
    // physics -> core, math, sim
    physics_mod.addImport("core", core_mod);
    physics_mod.addImport("math", math_mod);
    physics_mod.addImport("sim", sim_mod);
    // sim -> core, math (no physics needed unless referenced)
    sim_mod.addImport("core", core_mod);
    sim_mod.addImport("math", math_mod);
    // planet -> math
    planet_mod.addImport("math", math_mod);

    // wire imports so @import("core") etc. works in ffi
    root_mod.addImport("core", core_mod);
    root_mod.addImport("math", math_mod);
    root_mod.addImport("physics", physics_mod);
    root_mod.addImport("sim", sim_mod);
    root_mod.addImport("planet", planet_mod);

    // keep exports visible for wasm
    wasm_exe.entry = .disabled;
    wasm_exe.rdynamic = true;

    b.installArtifact(wasm_exe);
}

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
    // sim -> core, math
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

    // ----------------------------
    // Unit tests (native host)
    // ----------------------------
    const host_target = b.resolveTargetQuery(.{});
    const math_host = b.createModule(.{ .root_source_file = b.path("src/math/index.zig"), .target = host_target, .optimize = optimize });
    const core_host = b.createModule(.{ .root_source_file = b.path("src/core/types.zig"), .target = host_target, .optimize = optimize });
    const physics_host = b.createModule(.{ .root_source_file = b.path("src/physics/nbody.zig"), .target = host_target, .optimize = optimize });
    const sim_host = b.createModule(.{ .root_source_file = b.path("src/sim/index.zig"), .target = host_target, .optimize = optimize });
    const planet_host = b.createModule(.{ .root_source_file = b.path("src/planet/index.zig"), .target = host_target, .optimize = optimize });

    core_host.addImport("math", math_host);
    physics_host.addImport("core", core_host);
    physics_host.addImport("math", math_host);
    physics_host.addImport("sim", sim_host);
    sim_host.addImport("core", core_host);
    sim_host.addImport("math", math_host);
    planet_host.addImport("math", math_host);

    const math_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/math/index.zig"),
            .target = host_target,
            .optimize = optimize,
        }),
    });

    const grid_test_module = b.createModule(.{
        .root_source_file = b.path("src/planet/grid_test.zig"),
        .target = host_target,
        .optimize = optimize,
    });
    grid_test_module.addImport("planet", planet_host);
    const grid_tests = b.addTest(.{ .root_module = grid_test_module });

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&math_tests.step);
    test_step.dependOn(&grid_tests.step);

    // ----------------------------
    // Tool: Grid Generator (native host)
    // ----------------------------
    const gridgen_mod = b.createModule(.{
        .root_source_file = b.path("src/tools/gridgen.zig"),
        .target = host_target,
        .optimize = optimize,
    });
    const gridgen_exe = b.addExecutable(.{
        .name = "gridgen",
        .root_module = gridgen_mod,
    });
    // Reuse the host-native planet module which includes math
    gridgen_mod.addImport("planet", planet_host);

    b.installArtifact(gridgen_exe);

    const run_gridgen = b.addRunArtifact(gridgen_exe);
    const gridgen_step = b.step("gridgen", "Generate prebaked grid");
    gridgen_step.dependOn(&run_gridgen.step);
}

const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});

    const lib = b.addLibrary(.{
        .name = "simulation-service",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/lib.zig"),
            .target = target,
            .optimize = optimize,
        }),
        .linkage = .dynamic,
    });

    // Dependency: game-sim (for shader path if needed)
    const game_sim_dep = b.dependency("game_sim", .{});
    _ = game_sim_dep; // currently unused; shader is loaded at runtime

    // Shared C import module for wgpu-native
    const wgpu_module = b.createModule(.{ .root_source_file = b.path("src/wgpu_native.zig") });
    wgpu_module.addIncludePath(b.path("deps/wgpu-native/include"));
    lib.root_module.addImport("wgpu_native", wgpu_module);

    // Helper module depends on wgpu_native
    const helpers_module = b.createModule(.{ .root_source_file = b.path("src/wgpu_helpers.zig") });
    helpers_module.addImport("wgpu_native", wgpu_module);
    lib.root_module.addImport("wgpu_helpers", helpers_module);

    // Link wgpu-native from server-local deps
    lib.addIncludePath(b.path("deps/wgpu-native/include"));
    lib.addLibraryPath(b.path("deps/wgpu-native/lib"));
    lib.addRPath(b.path("deps/wgpu-native/lib"));
    lib.linkSystemLibrary("wgpu_native");

    b.installArtifact(lib);
}

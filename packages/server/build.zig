const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});

    const exe = b.addExecutable(.{
        .name = "native-host",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/native_host.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    // Shared C import module for wgpu-native
    const wgpu_module = b.createModule(.{ .root_source_file = b.path("src/wgpu_native.zig") });
    wgpu_module.addIncludePath(b.path("deps/wgpu-native/include"));
    exe.root_module.addImport("wgpu_native", wgpu_module);

    // Helper module depends on wgpu_native
    const helpers_module = b.createModule(.{ .root_source_file = b.path("src/wgpu_helpers.zig") });
    helpers_module.addImport("wgpu_native", wgpu_module);
    exe.root_module.addImport("wgpu_helpers", helpers_module);

    // Link wgpu-native from server-local deps
    exe.addIncludePath(b.path("deps/wgpu-native/include"));
    exe.addLibraryPath(b.path("deps/wgpu-native/lib"));
    exe.addRPath(b.path("deps/wgpu-native/lib"));
    exe.linkSystemLibrary("wgpu_native");

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    const run_step = b.step("run", "Run the native simulation host");
    run_step.dependOn(&run_cmd.step);
}

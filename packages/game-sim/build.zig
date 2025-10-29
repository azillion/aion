const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});

    //============================================================================
    // Configuration (wgpu-native via env vars)
    //============================================================================
    const wgpu_include_dir = std.process.getEnvVarOwned(b.allocator, "WGPU_NATIVE_INCLUDE_DIR") catch null;
    const wgpu_lib_dir = std.process.getEnvVarOwned(b.allocator, "WGPU_NATIVE_LIB_DIR") catch null;
    const wgpu_lib_name = std.process.getEnvVarOwned(b.allocator, "WGPU_NATIVE_LIB_NAME") catch "wgpu_native";

    var enable_headers: bool = false;
    if (wgpu_include_dir) |dir| {
        if (std.fs.path.join(b.allocator, &.{ dir, "webgpu", "wgpu.h" }) catch null) |p| {
            if (std.fs.cwd().openFile(p, .{})) |f| {
                f.close();
                enable_headers = true;
            } else |_| {}
            b.allocator.free(p);
        }
        if (!enable_headers) {
            if (std.fs.path.join(b.allocator, &.{ dir, "wgpu.h" }) catch null) |p2| {
                if (std.fs.cwd().openFile(p2, .{})) |f2| {
                    f2.close();
                    enable_headers = true;
                } else |_| {}
                b.allocator.free(p2);
            }
        }
    }

    const build_opts = b.addOptions();
    build_opts.addOption(bool, "enable_wgpu_headers", enable_headers);

    //============================================================================
    // Wasm Build -> Static Library
    //============================================================================
    const wasm_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    });

    const wasm_lib = b.addLibrary(.{
        .name = "game-sim",
        .linkage = .static,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = wasm_target,
            .optimize = optimize,
        }),
    });
    // Pass options to the module (0.15 API)
    wasm_lib.root_module.addOptions("build_options", build_opts);
    // Export the 'add' symbol for JavaScript to call

    if (enable_headers) if (wgpu_include_dir) |dir| wasm_lib.root_module.addIncludePath(b.path(dir));

    b.installArtifact(wasm_lib);

    //============================================================================
    // Native Build -> Shared Library
    //============================================================================
    const native_target = b.standardTargetOptions(.{});

    const native_lib = b.addLibrary(.{
        .name = "game-sim-native",
        .linkage = .dynamic,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = native_target,
            .optimize = optimize,
        }),
    });
    // Pass options to the module (0.15 API)
    native_lib.root_module.addOptions("build_options", build_opts);

    if (enable_headers) if (wgpu_include_dir) |dir| native_lib.root_module.addIncludePath(b.path(dir));

    if (wgpu_lib_dir) |libdir| {
        native_lib.addLibraryPath(b.path(libdir));
        native_lib.linkSystemLibrary(wgpu_lib_name);
    }

    b.installArtifact(native_lib);
}

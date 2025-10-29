const std = @import("std");
const wgpu = @import("wgpu_native").wgpu;

const helpers = @import("wgpu_helpers.zig");

inline fn sv(s: []const u8) wgpu.WGPUStringView {
    return .{ .data = @as([*c]const u8, @ptrCast(s.ptr)), .length = s.len };
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();

    std.debug.print("Initializing native WGPU device...\n", .{});

    const device = helpers.getHeadlessDevice() catch |err| {
        std.debug.print("Failed to initialize device: {s}\n", .{@errorName(err)});
        return err;
    };
    defer wgpu.wgpuDeviceRelease(device);

    const queue = wgpu.wgpuDeviceGetQueue(device);

    std.debug.print("Native device created. Setting up simulation...\n", .{});

    // Create ping-pong textures
    var grid_size: wgpu.WGPUExtent3D = .{ .width = 256, .height = 256, .depthOrArrayLayers = 1 };

    var tex_desc: wgpu.WGPUTextureDescriptor = .{};
    tex_desc.nextInChain = null;
    tex_desc.label = sv("GOL State Texture (Native)");
    tex_desc.usage = wgpu.WGPUTextureUsage_TextureBinding | wgpu.WGPUTextureUsage_StorageBinding | wgpu.WGPUTextureUsage_CopyDst;
    tex_desc.dimension = wgpu.WGPUTextureDimension_2D;
    tex_desc.size = grid_size;
    tex_desc.format = wgpu.WGPUTextureFormat_R32Uint;
    tex_desc.mipLevelCount = 1;
    tex_desc.sampleCount = 1;
    tex_desc.viewFormatCount = 0;
    tex_desc.viewFormats = null;

    const state_a = wgpu.wgpuDeviceCreateTexture(device, &tex_desc);
    defer wgpu.wgpuTextureRelease(state_a);
    const state_b = wgpu.wgpuDeviceCreateTexture(device, &tex_desc);
    defer wgpu.wgpuTextureRelease(state_b);

    const state_a_view = wgpu.wgpuTextureCreateView(state_a, null);
    defer wgpu.wgpuTextureViewRelease(state_a_view);
    const state_b_view = wgpu.wgpuTextureCreateView(state_b, null);
    defer wgpu.wgpuTextureViewRelease(state_b_view);

    // Seed random initial data into state_a
    var prng = std.Random.DefaultPrng.init(0xdeadbeef);
    const random = prng.random();
    var init_data: [256 * 256]u32 = undefined;
    var ii: usize = 0;
    while (ii < init_data.len) : (ii += 1) {
        init_data[ii] = if ((random.int(u32) & 3) == 0) 1 else 0;
    }
    const data_bytes = std.mem.sliceAsBytes(init_data[0..]);

    var dst: wgpu.WGPUTexelCopyTextureInfo = .{};
    dst.texture = state_a;
    dst.mipLevel = 0;
    dst.origin = .{ .x = 0, .y = 0, .z = 0 };
    dst.aspect = wgpu.WGPUTextureAspect_All;

    var layout: wgpu.WGPUTexelCopyBufferLayout = .{};
    layout.offset = 0;
    layout.bytesPerRow = 256 * 4;
    layout.rowsPerImage = 256;

    wgpu.wgpuQueueWriteTexture(queue, &dst, data_bytes.ptr, data_bytes.len, &layout, &grid_size);

    // Load compute shader WGSL at runtime from sibling package
    var shader_file = try std.fs.cwd().openFile("../game-sim/src/gol.wgsl", .{});
    defer shader_file.close();
    const source = try shader_file.readToEndAlloc(std.heap.page_allocator, 1024 * 1024);
    defer std.heap.page_allocator.free(source);
    var wgsl_src: wgpu.WGPUShaderSourceWGSL = .{};
    wgsl_src.chain.next = null;
    wgsl_src.chain.sType = wgpu.WGPUSType_ShaderSourceWGSL;
    wgsl_src.code = sv(source);

    var shader_desc: wgpu.WGPUShaderModuleDescriptor = .{};
    shader_desc.nextInChain = @as(*const wgpu.WGPUChainedStruct, @ptrCast(&wgsl_src.chain));
    shader_desc.label = sv("GOL Shader (Native)");

    const shader_module = wgpu.wgpuDeviceCreateShaderModule(device, &shader_desc);
    defer wgpu.wgpuShaderModuleRelease(shader_module);

    var bgl_entries: [2]wgpu.WGPUBindGroupLayoutEntry = .{ .{}, .{} };
    bgl_entries[0].nextInChain = null;
    bgl_entries[0].binding = 0;
    bgl_entries[0].visibility = wgpu.WGPUShaderStage_Compute;
    bgl_entries[0].storageTexture.access = wgpu.WGPUStorageTextureAccess_ReadOnly;
    bgl_entries[0].storageTexture.format = wgpu.WGPUTextureFormat_R32Uint;
    bgl_entries[0].storageTexture.viewDimension = wgpu.WGPUTextureViewDimension_2D;

    bgl_entries[1].nextInChain = null;
    bgl_entries[1].binding = 1;
    bgl_entries[1].visibility = wgpu.WGPUShaderStage_Compute;
    bgl_entries[1].storageTexture.access = wgpu.WGPUStorageTextureAccess_WriteOnly;
    bgl_entries[1].storageTexture.format = wgpu.WGPUTextureFormat_R32Uint;
    bgl_entries[1].storageTexture.viewDimension = wgpu.WGPUTextureViewDimension_2D;

    var bgl_desc: wgpu.WGPUBindGroupLayoutDescriptor = .{};
    bgl_desc.nextInChain = null;
    bgl_desc.label = sv("gol_bgl");
    bgl_desc.entryCount = bgl_entries.len;
    bgl_desc.entries = &bgl_entries;

    const bind_group_layout = wgpu.wgpuDeviceCreateBindGroupLayout(device, &bgl_desc);
    defer wgpu.wgpuBindGroupLayoutRelease(bind_group_layout);

    var pl_desc: wgpu.WGPUPipelineLayoutDescriptor = .{};
    pl_desc.nextInChain = null;
    pl_desc.label = sv("gol_pipeline_layout");
    pl_desc.bindGroupLayoutCount = 1;
    pl_desc.bindGroupLayouts = &bind_group_layout;
    const pipeline_layout = wgpu.wgpuDeviceCreatePipelineLayout(device, &pl_desc);
    defer wgpu.wgpuPipelineLayoutRelease(pipeline_layout);

    var stage: wgpu.WGPUProgrammableStageDescriptor = .{};
    stage.nextInChain = null;
    stage.module = shader_module;
    stage.entryPoint = sv("main");

    var cp_desc: wgpu.WGPUComputePipelineDescriptor = .{};
    cp_desc.nextInChain = null;
    cp_desc.label = sv("gol_compute_pipeline");
    cp_desc.layout = pipeline_layout;
    cp_desc.compute = stage;
    const pipeline = wgpu.wgpuDeviceCreateComputePipeline(device, &cp_desc);
    defer wgpu.wgpuComputePipelineRelease(pipeline);

    var entries_abw: [2]wgpu.WGPUBindGroupEntry = .{ .{}, .{} };
    entries_abw[0].nextInChain = null;
    entries_abw[0].binding = 0;
    entries_abw[0].textureView = state_a_view;
    entries_abw[1].nextInChain = null;
    entries_abw[1].binding = 1;
    entries_abw[1].textureView = state_b_view;

    var bg_desc_abw: wgpu.WGPUBindGroupDescriptor = .{};
    bg_desc_abw.nextInChain = null;
    bg_desc_abw.label = sv("bind_group_a_reads_b_writes");
    bg_desc_abw.layout = bind_group_layout;
    bg_desc_abw.entryCount = entries_abw.len;
    bg_desc_abw.entries = &entries_abw;
    const bind_group_a_reads_b_writes = wgpu.wgpuDeviceCreateBindGroup(device, &bg_desc_abw);
    defer wgpu.wgpuBindGroupRelease(bind_group_a_reads_b_writes);

    var entries_baw: [2]wgpu.WGPUBindGroupEntry = .{ .{}, .{} };
    entries_baw[0].nextInChain = null;
    entries_baw[0].binding = 0;
    entries_baw[0].textureView = state_b_view;
    entries_baw[1].nextInChain = null;
    entries_baw[1].binding = 1;
    entries_baw[1].textureView = state_a_view;

    var bg_desc_baw: wgpu.WGPUBindGroupDescriptor = .{};
    bg_desc_baw.nextInChain = null;
    bg_desc_baw.label = sv("bind_group_b_reads_a_writes");
    bg_desc_baw.layout = bind_group_layout;
    bg_desc_baw.entryCount = entries_baw.len;
    bg_desc_baw.entries = &entries_baw;
    const bind_group_b_reads_a_writes = wgpu.wgpuDeviceCreateBindGroup(device, &bg_desc_baw);
    defer wgpu.wgpuBindGroupRelease(bind_group_b_reads_a_writes);

    std.debug.print("Simulation setup complete. Ticking 100 times...\n", .{});
    var i: usize = 0;
    var is_state_a_current: bool = true;
    while (i < 100) : (i += 1) {
        var enc_desc: wgpu.WGPUCommandEncoderDescriptor = .{};
        enc_desc.nextInChain = null;
        enc_desc.label = sv("gol_cmd_encoder");
        const encoder = wgpu.wgpuDeviceCreateCommandEncoder(device, &enc_desc);

        var pass_desc: wgpu.WGPUComputePassDescriptor = .{};
        pass_desc.nextInChain = null;
        pass_desc.label = sv("gol_compute_pass");
        pass_desc.timestampWrites = null;
        const pass = wgpu.wgpuCommandEncoderBeginComputePass(encoder, &pass_desc);
        wgpu.wgpuComputePassEncoderSetPipeline(pass, pipeline);
        if (is_state_a_current) {
            wgpu.wgpuComputePassEncoderSetBindGroup(pass, 0, bind_group_a_reads_b_writes, 0, null);
        } else {
            wgpu.wgpuComputePassEncoderSetBindGroup(pass, 0, bind_group_b_reads_a_writes, 0, null);
        }
        // 256x256 grid with 8x8 WG size => 32x32 workgroups
        wgpu.wgpuComputePassEncoderDispatchWorkgroups(pass, 32, 32, 1);
        wgpu.wgpuComputePassEncoderEnd(pass);

        const cmd = wgpu.wgpuCommandEncoderFinish(encoder, null);
        wgpu.wgpuQueueSubmit(queue, 1, &cmd);

        _ = wgpu.wgpuDevicePoll(device, 1, null);

        is_state_a_current = !is_state_a_current;
        wgpu.wgpuComputePassEncoderRelease(pass);
        wgpu.wgpuCommandEncoderRelease(encoder);
        wgpu.wgpuCommandBufferRelease(cmd);
    }

    std.debug.print("Native simulation finished successfully.\n", .{});
}

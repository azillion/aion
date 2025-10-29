const std = @import("std");
const build_options = @import("build_options");
const wgpu = blk: {
    if (build_options.enable_wgpu_headers) {
        break :blk @cImport(@cInclude("webgpu/webgpu.h"));
    } else {
        break :blk struct {};
    }
};

inline fn sv(s: []const u8) wgpu.struct_WGPUStringView {
    return .{ .data = @as([*c]const u8, @ptrCast(s.ptr)), .length = s.len };
}

const SimulatorState = struct {
    device: wgpu.WGPUDevice,
    queue: wgpu.WGPUQueue,
    allocator: std.mem.Allocator,

    state_a: wgpu.WGPUTexture,
    state_b: wgpu.WGPUTexture,
    state_a_view: wgpu.WGPUTextureView,
    state_b_view: wgpu.WGPUTextureView,

    is_state_a_current: bool,
    grid_size: wgpu.WGPUExtent3D,

    bind_group_layout: wgpu.WGPUBindGroupLayout,
    pipeline: wgpu.WGPUComputePipeline,

    bind_group_a_reads_b_writes: wgpu.WGPUBindGroup,
    bind_group_b_reads_a_writes: wgpu.WGPUBindGroup,
};

export fn create_simulator(device: wgpu.WGPUDevice, queue: wgpu.WGPUQueue) callconv(.c) ?*SimulatorState {
    const builtin = @import("builtin");
    const allocator: std.mem.Allocator = if (builtin.target.cpu.arch == .wasm32) std.heap.wasm_allocator else std.heap.page_allocator;

    const state = allocator.create(SimulatorState) catch return null;
    state.* = .{
        .device = device,
        .queue = queue,
        .allocator = allocator,
        .state_a = null,
        .state_b = null,
        .state_a_view = null,
        .state_b_view = null,
        .is_state_a_current = true,
        .grid_size = .{ .width = 256, .height = 256, .depthOrArrayLayers = 1 },
        .bind_group_layout = null,
        .pipeline = null,
        .bind_group_a_reads_b_writes = null,
        .bind_group_b_reads_a_writes = null,
    };

    // Create ping-pong textures
    var tex_desc: wgpu.WGPUTextureDescriptor = .{};
    tex_desc.nextInChain = null;
    tex_desc.label = sv("gol_state_texture");
    tex_desc.usage = wgpu.WGPUTextureUsage_TextureBinding | wgpu.WGPUTextureUsage_StorageBinding | wgpu.WGPUTextureUsage_CopyDst;
    tex_desc.dimension = wgpu.WGPUTextureDimension_2D;
    tex_desc.size = state.grid_size;
    tex_desc.format = wgpu.WGPUTextureFormat_R32Uint;
    tex_desc.mipLevelCount = 1;
    tex_desc.sampleCount = 1;
    tex_desc.viewFormatCount = 0;
    tex_desc.viewFormats = null;

    state.state_a = wgpu.wgpuDeviceCreateTexture(state.device, &tex_desc);
    state.state_b = wgpu.wgpuDeviceCreateTexture(state.device, &tex_desc);

    state.state_a_view = wgpu.wgpuTextureCreateView(state.state_a, null);
    state.state_b_view = wgpu.wgpuTextureCreateView(state.state_b, null);

    // Seed initial random data into state_a
    var prng = std.Random.DefaultPrng.init(0xdeadbeef);
    const random = prng.random();
    var init_data: [256 * 256]u32 = undefined;
    var i: usize = 0;
    while (i < init_data.len) : (i += 1) {
        // 25% chance alive
        const alive: u32 = if ((random.int(u32) & 3) == 0) 1 else 0;
        init_data[i] = alive;
    }

    const data_slice: []const u32 = init_data[0..];
    const data_bytes: []const u8 = std.mem.sliceAsBytes(data_slice);

    var dst: wgpu.WGPUTexelCopyTextureInfo = .{};
    dst.texture = state.state_a;
    dst.mipLevel = 0;
    dst.origin = .{ .x = 0, .y = 0, .z = 0 };
    dst.aspect = wgpu.WGPUTextureAspect_All;

    var layout: wgpu.WGPUTexelCopyBufferLayout = .{};
    layout.offset = 0;
    layout.bytesPerRow = 256 * 4;
    layout.rowsPerImage = 256;

    wgpu.wgpuQueueWriteTexture(
        state.queue,
        &dst,
        data_bytes.ptr,
        data_bytes.len,
        &layout,
        &state.grid_size,
    );

    // Create compute pipeline
    const source = @embedFile("gol.wgsl");

    var wgsl_src: wgpu.WGPUShaderSourceWGSL = .{};
    wgsl_src.chain.next = null;
    wgsl_src.chain.sType = wgpu.WGPUSType_ShaderSourceWGSL;
    wgsl_src.code = sv(source);

    var shader_desc: wgpu.WGPUShaderModuleDescriptor = .{};
    shader_desc.nextInChain = @as(*const wgpu.WGPUChainedStruct, @ptrCast(&wgsl_src.chain));
    shader_desc.label = sv("gol_wgsl");

    const shader_module = wgpu.wgpuDeviceCreateShaderModule(state.device, &shader_desc);

    var bgl_entries: [2]wgpu.WGPUBindGroupLayoutEntry = .{ .{}, .{} };

    // Binding 0: read-only storage texture (input)
    bgl_entries[0].nextInChain = null;
    bgl_entries[0].binding = 0;
    bgl_entries[0].visibility = wgpu.WGPUShaderStage_Compute;
    bgl_entries[0].storageTexture.access = wgpu.WGPUStorageTextureAccess_ReadOnly;
    bgl_entries[0].storageTexture.format = wgpu.WGPUTextureFormat_R32Uint;
    bgl_entries[0].storageTexture.viewDimension = wgpu.WGPUTextureViewDimension_2D;

    // Binding 1: write-only storage texture (output)
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

    state.bind_group_layout = wgpu.wgpuDeviceCreateBindGroupLayout(state.device, &bgl_desc);

    var pl_desc: wgpu.WGPUPipelineLayoutDescriptor = .{};
    pl_desc.nextInChain = null;
    pl_desc.label = sv("gol_pipeline_layout");
    pl_desc.bindGroupLayoutCount = 1;
    pl_desc.bindGroupLayouts = &state.bind_group_layout;

    const pipeline_layout = wgpu.wgpuDeviceCreatePipelineLayout(state.device, &pl_desc);

    var stage: wgpu.WGPUProgrammableStageDescriptor = .{};
    stage.nextInChain = null;
    stage.module = shader_module;
    stage.entryPoint = sv("main");

    var cp_desc: wgpu.WGPUComputePipelineDescriptor = .{};
    cp_desc.nextInChain = null;
    cp_desc.label = sv("gol_compute_pipeline");
    cp_desc.layout = pipeline_layout;
    cp_desc.compute = stage;

    state.pipeline = wgpu.wgpuDeviceCreateComputePipeline(state.device, &cp_desc);

    // Shader module and pipeline layout no longer needed after pipeline creation
    wgpu.wgpuShaderModuleRelease(shader_module);
    wgpu.wgpuPipelineLayoutRelease(pipeline_layout);

    // Create bind groups for ping-pong
    var entries_abw: [2]wgpu.WGPUBindGroupEntry = .{ .{}, .{} };
    entries_abw[0].nextInChain = null;
    entries_abw[0].binding = 0;
    entries_abw[0].textureView = state.state_a_view;
    entries_abw[1].nextInChain = null;
    entries_abw[1].binding = 1;
    entries_abw[1].textureView = state.state_b_view;

    var bg_desc_abw: wgpu.WGPUBindGroupDescriptor = .{};
    bg_desc_abw.nextInChain = null;
    bg_desc_abw.label = sv("bind_group_a_reads_b_writes");
    bg_desc_abw.layout = state.bind_group_layout;
    bg_desc_abw.entryCount = entries_abw.len;
    bg_desc_abw.entries = &entries_abw;

    state.bind_group_a_reads_b_writes = wgpu.wgpuDeviceCreateBindGroup(state.device, &bg_desc_abw);

    var entries_baw: [2]wgpu.WGPUBindGroupEntry = .{ .{}, .{} };
    entries_baw[0].nextInChain = null;
    entries_baw[0].binding = 0;
    entries_baw[0].textureView = state.state_b_view;
    entries_baw[1].nextInChain = null;
    entries_baw[1].binding = 1;
    entries_baw[1].textureView = state.state_a_view;

    var bg_desc_baw: wgpu.WGPUBindGroupDescriptor = .{};
    bg_desc_baw.nextInChain = null;
    bg_desc_baw.label = sv("bind_group_b_reads_a_writes");
    bg_desc_baw.layout = state.bind_group_layout;
    bg_desc_baw.entryCount = entries_baw.len;
    bg_desc_baw.entries = &entries_baw;

    state.bind_group_b_reads_a_writes = wgpu.wgpuDeviceCreateBindGroup(state.device, &bg_desc_baw);

    return state;
}

export fn tick_simulator(state: *SimulatorState) callconv(.c) void {
    var enc_desc: wgpu.WGPUCommandEncoderDescriptor = .{};
    enc_desc.nextInChain = null;
    enc_desc.label = sv("gol_cmd_encoder");

    const encoder = wgpu.wgpuDeviceCreateCommandEncoder(state.device, &enc_desc);

    var pass_desc: wgpu.WGPUComputePassDescriptor = .{};
    pass_desc.nextInChain = null;
    pass_desc.label = sv("gol_compute_pass");
    pass_desc.timestampWrites = null;

    const pass = wgpu.wgpuCommandEncoderBeginComputePass(encoder, &pass_desc);
    wgpu.wgpuComputePassEncoderSetPipeline(pass, state.pipeline);

    if (state.is_state_a_current) {
        wgpu.wgpuComputePassEncoderSetBindGroup(pass, 0, state.bind_group_a_reads_b_writes, 0, null);
    } else {
        wgpu.wgpuComputePassEncoderSetBindGroup(pass, 0, state.bind_group_b_reads_a_writes, 0, null);
    }

    // 256x256 grid with 8x8 workgroup size => 32x32 workgroups
    wgpu.wgpuComputePassEncoderDispatchWorkgroups(pass, 32, 32, 1);
    wgpu.wgpuComputePassEncoderEnd(pass);

    const cmd = wgpu.wgpuCommandEncoderFinish(encoder, null);
    wgpu.wgpuQueueSubmit(state.queue, 1, &cmd);

    // Optionally release temporary objects
    wgpu.wgpuComputePassEncoderRelease(pass);
    wgpu.wgpuCommandEncoderRelease(encoder);
    wgpu.wgpuCommandBufferRelease(cmd);

    state.is_state_a_current = !state.is_state_a_current;
}

export fn get_output_texture_view(state: *SimulatorState) callconv(.c) wgpu.WGPUTextureView {
    if (state.is_state_a_current) {
        return state.state_b_view;
    } else {
        return state.state_a_view;
    }
}

export fn destroy_simulator(state: *SimulatorState) callconv(.c) void {
    // Release GPU resources first
    if (state.state_a_view != null) wgpu.wgpuTextureViewRelease(state.state_a_view);
    if (state.state_b_view != null) wgpu.wgpuTextureViewRelease(state.state_b_view);
    if (state.state_a != null) wgpu.wgpuTextureRelease(state.state_a);
    if (state.state_b != null) wgpu.wgpuTextureRelease(state.state_b);

    if (state.bind_group_a_reads_b_writes != null) wgpu.wgpuBindGroupRelease(state.bind_group_a_reads_b_writes);
    if (state.bind_group_b_reads_a_writes != null) wgpu.wgpuBindGroupRelease(state.bind_group_b_reads_a_writes);

    if (state.pipeline != null) wgpu.wgpuComputePipelineRelease(state.pipeline);
    if (state.bind_group_layout != null) wgpu.wgpuBindGroupLayoutRelease(state.bind_group_layout);

    const allocator = state.allocator;
    allocator.destroy(state);
}

const std = @import("std");
const wgpu = @import("wgpu_native").wgpu;
const helpers = @import("wgpu_helpers.zig");

inline fn sv(s: []const u8) wgpu.WGPUStringView {
    return .{ .data = @as([*c]const u8, @ptrCast(s.ptr)), .length = s.len };
}

pub const SimState = struct {
    device: wgpu.WGPUDevice,
    queue: wgpu.WGPUQueue,
    pipeline: wgpu.WGPUComputePipeline,
    state_a: wgpu.WGPUTexture,
    state_b: wgpu.WGPUTexture,
    state_a_view: wgpu.WGPUTextureView,
    state_b_view: wgpu.WGPUTextureView,
    bg_a_reads_b: wgpu.WGPUBindGroup,
    bg_b_reads_a: wgpu.WGPUBindGroup,
    bind_group_layout: wgpu.WGPUBindGroupLayout,
    readback_buffer: wgpu.WGPUBuffer,
    is_state_a_current: bool,
    grid_size: wgpu.WGPUExtent3D,
    buffer_size: usize,
};

// Map callback support
const MapInfo = struct {
    finished: bool = false,
    status: wgpu.WGPUMapAsyncStatus = wgpu.WGPUMapAsyncStatus_Success,
};

fn mapCallback(status: wgpu.WGPUMapAsyncStatus, _: wgpu.WGPUStringView, userdata1: ?*anyopaque, _: ?*anyopaque) callconv(.c) void {
    if (userdata1) |ud| {
        var info: *MapInfo = @ptrCast(@alignCast(ud));
        info.status = status;
        info.finished = true;
    }
}

export fn create_simulator() ?*SimState {
    const device = helpers.getHeadlessDevice() catch return null;
    const queue = wgpu.wgpuDeviceGetQueue(device);

    var grid_size: wgpu.WGPUExtent3D = .{ .width = 2048, .height = 2048, .depthOrArrayLayers = 1 };
    const buffer_size: usize = grid_size.width * grid_size.height * 4;

    // Create ping-pong textures
    var tex_desc: wgpu.WGPUTextureDescriptor = .{};
    tex_desc.nextInChain = null;
    tex_desc.label = sv("GOL State Texture");
    tex_desc.usage = wgpu.WGPUTextureUsage_TextureBinding | wgpu.WGPUTextureUsage_StorageBinding | wgpu.WGPUTextureUsage_CopyDst | wgpu.WGPUTextureUsage_CopySrc;
    tex_desc.dimension = wgpu.WGPUTextureDimension_2D;
    tex_desc.size = grid_size;
    tex_desc.format = wgpu.WGPUTextureFormat_R32Uint;
    tex_desc.mipLevelCount = 1;
    tex_desc.sampleCount = 1;
    tex_desc.viewFormatCount = 0;
    tex_desc.viewFormats = null;

    const state_a = wgpu.wgpuDeviceCreateTexture(device, &tex_desc);
    const state_b = wgpu.wgpuDeviceCreateTexture(device, &tex_desc);

    const state_a_view = wgpu.wgpuTextureCreateView(state_a, null);
    const state_b_view = wgpu.wgpuTextureCreateView(state_b, null);

    // Seed random initial data into state_a
    var prng = std.Random.DefaultPrng.init(0xdeadbeef);
    const random = prng.random();
    const cell_count: usize = @as(usize, grid_size.width) * @as(usize, grid_size.height);
    const init_data = std.heap.page_allocator.alloc(u32, cell_count) catch return null;
    defer std.heap.page_allocator.free(init_data);
    var ii: usize = 0;
    while (ii < cell_count) : (ii += 1) {
        init_data[ii] = if ((random.int(u32) & 3) == 0) 1 else 0;
    }
    const data_bytes = std.mem.sliceAsBytes(init_data);

    var dst: wgpu.WGPUTexelCopyTextureInfo = .{};
    dst.texture = state_a;
    dst.mipLevel = 0;
    dst.origin = .{ .x = 0, .y = 0, .z = 0 };
    dst.aspect = wgpu.WGPUTextureAspect_All;

    var layout: wgpu.WGPUTexelCopyBufferLayout = .{};
    layout.offset = 0;
    layout.bytesPerRow = grid_size.width * 4;
    layout.rowsPerImage = grid_size.height;

    wgpu.wgpuQueueWriteTexture(queue, &dst, data_bytes.ptr, data_bytes.len, &layout, &grid_size);

    // Load compute shader WGSL from sibling package at runtime
    var shader_file = std.fs.cwd().openFile("../game-sim/src/gol.wgsl", .{}) catch return null;
    defer shader_file.close();
    const source = shader_file.readToEndAlloc(std.heap.page_allocator, 1024 * 1024) catch return null;
    defer std.heap.page_allocator.free(source);
    var wgsl_src: wgpu.WGPUShaderSourceWGSL = .{};
    wgsl_src.chain.next = null;
    wgsl_src.chain.sType = wgpu.WGPUSType_ShaderSourceWGSL;
    wgsl_src.code = sv(source);

    var shader_desc: wgpu.WGPUShaderModuleDescriptor = .{};
    shader_desc.nextInChain = @as(*const wgpu.WGPUChainedStruct, @ptrCast(&wgsl_src.chain));
    shader_desc.label = sv("GOL Shader");
    const shader_module = wgpu.wgpuDeviceCreateShaderModule(device, &shader_desc);

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

    var pl_desc: wgpu.WGPUPipelineLayoutDescriptor = .{};
    pl_desc.nextInChain = null;
    pl_desc.label = sv("gol_pipeline_layout");
    pl_desc.bindGroupLayoutCount = 1;
    pl_desc.bindGroupLayouts = &bind_group_layout;
    const pipeline_layout = wgpu.wgpuDeviceCreatePipelineLayout(device, &pl_desc);

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

    // NOTE: The correct type is WGPUBindGroupEntry; use exact type
    // We split the declaration and assignment to avoid typos
    var entries_baw_fixed: [2]wgpu.WGPUBindGroupEntry = .{ .{}, .{} };
    entries_baw_fixed[0].nextInChain = null;
    entries_baw_fixed[0].binding = 0;
    entries_baw_fixed[0].textureView = state_b_view;
    entries_baw_fixed[1].nextInChain = null;
    entries_baw_fixed[1].binding = 1;
    entries_baw_fixed[1].textureView = state_a_view;

    var bg_desc_baw: wgpu.WGPUBindGroupDescriptor = .{};
    bg_desc_baw.nextInChain = null;
    bg_desc_baw.label = sv("bind_group_b_reads_a_writes");
    bg_desc_baw.layout = bind_group_layout;
    bg_desc_baw.entryCount = entries_baw_fixed.len;
    bg_desc_baw.entries = &entries_baw_fixed;
    const bind_group_b_reads_a_writes = wgpu.wgpuDeviceCreateBindGroup(device, &bg_desc_baw);

    const readback_buffer = wgpu.wgpuDeviceCreateBuffer(device, &.{
        .nextInChain = null,
        .label = sv("GOL Readback Buffer"),
        .usage = wgpu.WGPUBufferUsage_CopyDst | wgpu.WGPUBufferUsage_MapRead,
        .size = buffer_size,
        .mappedAtCreation = 0,
    });

    const state = std.heap.c_allocator.create(SimState) catch return null;
    state.* = .{
        .device = device,
        .queue = queue,
        .pipeline = pipeline,
        .state_a = state_a,
        .state_b = state_b,
        .state_a_view = state_a_view,
        .state_b_view = state_b_view,
        .bg_a_reads_b = bind_group_a_reads_b_writes,
        .bg_b_reads_a = bind_group_b_reads_a_writes,
        .bind_group_layout = bind_group_layout,
        .readback_buffer = readback_buffer,
        .is_state_a_current = true,
        .grid_size = grid_size,
        .buffer_size = buffer_size,
    };

    // Release temporary handles retained elsewhere
    wgpu.wgpuShaderModuleRelease(shader_module);
    wgpu.wgpuPipelineLayoutRelease(pipeline_layout);

    return state;
}

export fn tick_simulator(state: *SimState) void {
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
        wgpu.wgpuComputePassEncoderSetBindGroup(pass, 0, state.bg_a_reads_b, 0, null);
    } else {
        wgpu.wgpuComputePassEncoderSetBindGroup(pass, 0, state.bg_b_reads_a, 0, null);
    }
    // 256x256 grid with 8x8 WG size => 32x32 workgroups
    wgpu.wgpuComputePassEncoderDispatchWorkgroups(pass, state.grid_size.width / 8, state.grid_size.height / 8, 1);
    wgpu.wgpuComputePassEncoderEnd(pass);

    const cmd = wgpu.wgpuCommandEncoderFinish(encoder, null);
    wgpu.wgpuQueueSubmit(state.queue, 1, &cmd);

    state.is_state_a_current = !state.is_state_a_current;

    wgpu.wgpuComputePassEncoderRelease(pass);
    wgpu.wgpuCommandEncoderRelease(encoder);
    wgpu.wgpuCommandBufferRelease(cmd);
}

export fn get_simulation_data(state: *SimState) ?[*]u8 {
    const final_texture = if (state.is_state_a_current) state.state_a else state.state_b;
    var encoder_desc: wgpu.WGPUCommandEncoderDescriptor = .{};
    encoder_desc.nextInChain = null;
    encoder_desc.label = sv("readback_encoder");
    const encoder = wgpu.wgpuDeviceCreateCommandEncoder(state.device, &encoder_desc);

    var src_tex: wgpu.WGPUTexelCopyTextureInfo = .{};
    src_tex.texture = final_texture;
    src_tex.mipLevel = 0;
    src_tex.origin = .{ .x = 0, .y = 0, .z = 0 };
    src_tex.aspect = wgpu.WGPUTextureAspect_All;

    var dst_buf: wgpu.WGPUTexelCopyBufferInfo = .{};
    dst_buf.buffer = state.readback_buffer;
    dst_buf.layout = .{ .offset = 0, .bytesPerRow = state.grid_size.width * 4, .rowsPerImage = state.grid_size.height };

    wgpu.wgpuCommandEncoderCopyTextureToBuffer(encoder, &src_tex, &dst_buf, &state.grid_size);
    const cmd = wgpu.wgpuCommandEncoderFinish(encoder, null);
    wgpu.wgpuQueueSubmit(state.queue, 1, &cmd);

    var info = MapInfo{};
    var cbinfo: wgpu.WGPUBufferMapCallbackInfo = .{};
    cbinfo.nextInChain = null;
    cbinfo.mode = wgpu.WGPUCallbackMode_AllowSpontaneous;
    cbinfo.callback = mapCallback;
    cbinfo.userdata1 = &info;
    cbinfo.userdata2 = null;
    _ = wgpu.wgpuBufferMapAsync(state.readback_buffer, wgpu.WGPUMapMode_Read, 0, state.buffer_size, cbinfo);
    while (!info.finished) {
        _ = wgpu.wgpuDevicePoll(state.device, 1, null);
    }

    if (info.status != wgpu.WGPUMapAsyncStatus_Success) {
        std.debug.print("Zig: MAP FAILED! Returning null.\n", .{});
        wgpu.wgpuCommandBufferRelease(cmd);
        wgpu.wgpuCommandEncoderRelease(encoder);
        return null;
    }

    const raw_ptr = wgpu.wgpuBufferGetMappedRange(state.readback_buffer, 0, state.buffer_size);
    const addr: usize = @intFromPtr(raw_ptr);
    if (addr == 0) {
        wgpu.wgpuCommandBufferRelease(cmd);
        wgpu.wgpuCommandEncoderRelease(encoder);
        return null;
    }

    const as_bytes: [*]u8 = @ptrCast(raw_ptr);

    wgpu.wgpuCommandBufferRelease(cmd);
    wgpu.wgpuCommandEncoderRelease(encoder);
    return as_bytes;
}

export fn release_simulation_data(state: *SimState) void {
    wgpu.wgpuBufferUnmap(state.readback_buffer);
}

export fn destroy_simulator(state: *SimState) void {
    // Release GPU resources
    wgpu.wgpuBufferDestroy(state.readback_buffer);
    wgpu.wgpuBufferRelease(state.readback_buffer);
    wgpu.wgpuBindGroupRelease(state.bg_a_reads_b);
    wgpu.wgpuBindGroupRelease(state.bg_b_reads_a);
    wgpu.wgpuBindGroupLayoutRelease(state.bind_group_layout);
    wgpu.wgpuTextureViewRelease(state.state_a_view);
    wgpu.wgpuTextureViewRelease(state.state_b_view);
    wgpu.wgpuTextureRelease(state.state_a);
    wgpu.wgpuTextureRelease(state.state_b);
    wgpu.wgpuComputePipelineRelease(state.pipeline);
    wgpu.wgpuDeviceRelease(state.device);

    std.heap.c_allocator.destroy(state);
}

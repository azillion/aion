const std = @import("std");
const wgpu = @import("wgpu_native").wgpu;

inline fn sv(s: []const u8) wgpu.WGPUStringView {
    return .{ .data = @as([*c]const u8, @ptrCast(s.ptr)), .length = s.len };
}

var g_adapter: wgpu.WGPUAdapter = null;
var g_device: wgpu.WGPUDevice = null;

fn requestAdapterCallback(status: wgpu.WGPURequestAdapterStatus, adapter: wgpu.WGPUAdapter, message: wgpu.WGPUStringView, _: ?*anyopaque, _: ?*anyopaque) callconv(.c) void {
    switch (status) {
        wgpu.WGPURequestAdapterStatus_Success => {
            g_adapter = adapter;
        },
        else => {
            const slice = @as([*]const u8, @ptrCast(message.data))[0..message.length];
            std.debug.print("Failed to get adapter: {s}\n", .{slice});
            g_adapter = null;
        },
    }
}

fn requestDeviceCallback(status: wgpu.WGPURequestDeviceStatus, device: wgpu.WGPUDevice, message: wgpu.WGPUStringView, _: ?*anyopaque, _: ?*anyopaque) callconv(.c) void {
    switch (status) {
        wgpu.WGPURequestDeviceStatus_Success => {
            g_device = device;
        },
        else => {
            const slice = @as([*]const u8, @ptrCast(message.data))[0..message.length];
            std.debug.print("Failed to get device: {s}\n", .{slice});
            g_device = null;
        },
    }
}

pub const InitError = error{ AdapterUnavailable, DeviceUnavailable, Timeout };

pub fn getHeadlessDevice() InitError!wgpu.WGPUDevice {
    const instance = wgpu.wgpuCreateInstance(null);
    if (instance == null) return InitError.AdapterUnavailable;

    g_adapter = null;
    const adapter_cb: wgpu.WGPURequestAdapterCallbackInfo = .{
        .nextInChain = null,
        .mode = wgpu.WGPUCallbackMode_AllowSpontaneous,
        .callback = requestAdapterCallback,
        .userdata1 = null,
        .userdata2 = null,
    };
    _ = wgpu.wgpuInstanceRequestAdapter(instance, null, adapter_cb);

    var spin: usize = 0;
    while (g_adapter == null and spin < 10_000_000) : (spin += 1) {}
    if (g_adapter == null) return InitError.Timeout;
    const adapter = g_adapter;

    g_device = null;
    const device_cb: wgpu.WGPURequestDeviceCallbackInfo = .{
        .nextInChain = null,
        .mode = wgpu.WGPUCallbackMode_AllowSpontaneous,
        .callback = requestDeviceCallback,
        .userdata1 = null,
        .userdata2 = null,
    };
    _ = wgpu.wgpuAdapterRequestDevice(adapter, null, device_cb);

    spin = 0;
    while (g_device == null and spin < 10_000_000) : (spin += 1) {}
    if (g_device == null) return InitError.DeviceUnavailable;

    wgpu.wgpuAdapterRelease(adapter);
    wgpu.wgpuInstanceRelease(instance);
    return g_device;
}

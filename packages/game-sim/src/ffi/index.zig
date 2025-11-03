const std = @import("std");
const sim = @import("sim");
const types = @import("core");
const math = @import("math");
const physics = @import("physics");
const planet = @import("planet");

// Import JS logging function for error reporting
extern "env" fn log_error(ptr: [*]const u8, len: usize) void;

// C-compatible slice return for WASM FFI
pub const SliceU8 = extern struct { ptr: u32, len: u32 };

// Temporary schema for adding a body from JSON (no id)
const JsonBody = struct {
    name: []const u8,
    position: types.Vec3,
    velocity: types.Vec3,
    radius: f64,
    mass: f64,
    albedo: types.Vec3,
    emissive: ?types.Vec3 = null,
};

fn wasmAllocator() std.mem.Allocator {
    const builtin = @import("builtin");
    return if (builtin.target.cpu.arch == .wasm32) std.heap.wasm_allocator else std.heap.page_allocator;
}

// Input buffer owned by the module; host writes mask here via get_input_buffer
var input_buffer: types.InputState = .{ .keys_mask = 0 };

// Scratch area for writing [ptr,len] results back to the host.
var query_scratch: [2]u32 = .{ 0, 0 };

pub export fn get_query_scratch_ptr() callconv(.c) u32 {
    return @intCast(@intFromPtr(&query_scratch));
}

// General purpose scratch buffer for string/JSON inputs from host
var scratch_buffer: [4096]u8 = undefined;
pub export fn get_scratch_buffer_ptr() callconv(.c) u32 {
    return @intCast(@intFromPtr(&scratch_buffer));
}

pub export fn get_input_buffer() callconv(.c) u32 {
    return @intCast(@intFromPtr(&input_buffer));
}

pub export fn create_simulator(initial_state_ptr: u32, initial_state_len: u32) callconv(.c) ?*sim.Simulator {
    const allocator = wasmAllocator();
    const ptr: [*]const u8 = @ptrFromInt(initial_state_ptr);
    const input_slice = ptr[0..initial_state_len];
    var parsed = std.json.parseFromSlice(types.SystemState, allocator, input_slice, .{ .ignore_unknown_fields = true }) catch |err| {
        var buf: [256]u8 = undefined;
        const msg = std.fmt.bufPrint(&buf, "JSON parse error: {s}", .{@errorName(err)}) catch {
            const lit = "JSON parse error";
            log_error(lit.ptr, lit.len);
            return null;
        };
        log_error(msg.ptr, msg.len);
        return null;
    };
    // Take ownership of parsed memory: do NOT deinit parsed here; hand it to Simulator
    const owned_state = parsed.value;
    const simulator = sim.Simulator.create(allocator, owned_state) catch {
        parsed.deinit();
        return null;
    };
    // Intentionally do not parsed.deinit(); simulator owns the memory
    return simulator;
}

pub export fn tick_simulator(simulator: *sim.Simulator, dt_unscaled: f64) callconv(.c) void {
    physics.tickSimulator(simulator, dt_unscaled, &input_buffer);
}

pub export fn query_simulator_state(simulator: *sim.Simulator) callconv(.c) SliceU8 {
    const alloc = wasmAllocator();
    // Use ArrayListUnmanaged; provide allocator at usage sites
    var a = std.ArrayListUnmanaged(u8){};
    defer a.deinit(alloc);

    var old = a.writer(alloc);
    var ad = old.adaptToNewApi(&.{});
    const w: *std.Io.Writer = &ad.new_interface;

    const fmt = std.json.fmt(simulator.state, .{});
    fmt.format(w) catch return .{ .ptr = 0, .len = 0 };

    const buf = a.toOwnedSlice(alloc) catch return .{ .ptr = 0, .len = 0 };
    return .{ .ptr = @intCast(@intFromPtr(buf.ptr)), .len = @intCast(buf.len) };
}

// Out-parameter variant to avoid ABI differences across engines.
pub export fn query_simulator_state_out(simulator: *sim.Simulator, out_ptr: u32) callconv(.c) void {
    const alloc = wasmAllocator();

    var a = std.ArrayListUnmanaged(u8){};
    defer a.deinit(alloc);

    var old = a.writer(alloc);
    var ad = old.adaptToNewApi(&.{});
    const w: *std.Io.Writer = &ad.new_interface;

    const fmt = std.json.fmt(simulator.state, .{});
    fmt.format(w) catch {
        const out: [*]u32 = @ptrFromInt(out_ptr);
        out[0] = 0;
        out[1] = 0;
        return;
    };

    const buf = a.toOwnedSlice(alloc) catch {
        const out: [*]u32 = @ptrFromInt(out_ptr);
        out[0] = 0;
        out[1] = 0;
        return;
    };

    const out: [*]u32 = @ptrFromInt(out_ptr);
    out[0] = @intCast(@intFromPtr(buf.ptr));
    out[1] = @intCast(buf.len);
}

// --- Control and action functions ---
pub export fn set_time_scale(simulator: *sim.Simulator, scale: f64) callconv(.c) void {
    simulator.time_scale = if (scale < 0.0) 0.0 else scale;
}

pub export fn add_body(simulator: *sim.Simulator, body_json_ptr: u32, body_json_len: u32) callconv(.c) void {
    const allocator = simulator.allocator;
    const ptr: [*]const u8 = @ptrFromInt(body_json_ptr);
    const json_slice = ptr[0..body_json_len];

    var parsed = std.json.parseFromSlice(JsonBody, allocator, json_slice, .{ .ignore_unknown_fields = true }) catch {
        return;
    };
    defer parsed.deinit();

    const name_dup = allocator.dupe(u8, parsed.value.name) catch {
        return;
    };
    const id_str = std.fmt.allocPrint(allocator, "body-{d}", .{simulator.next_body_id}) catch {
        allocator.free(name_dup);
        return;
    };
    simulator.next_body_id += 1;

    const new_body = types.Body{
        .id = id_str,
        .name = name_dup,
        .position = parsed.value.position,
        .velocity = parsed.value.velocity,
        .radius = parsed.value.radius,
        .mass = parsed.value.mass,
        .albedo = parsed.value.albedo,
        .emissive = parsed.value.emissive,
        .terrain = null,
    };

    const old_len = simulator.state.bodies.len;
    const new_slice = allocator.realloc(simulator.state.bodies, old_len + 1) catch {
        allocator.free(id_str);
        allocator.free(name_dup);
        return;
    };
    new_slice[old_len] = new_body;
    simulator.state.bodies = new_slice;
}

pub export fn teleport_to_surface(simulator: *sim.Simulator, target_id_ptr: u32, target_id_len: u32) callconv(.c) void {
    const id_ptr: [*]const u8 = @ptrFromInt(target_id_ptr);
    const target_id = id_ptr[0..target_id_len];
    const player_id = "player-ship";

    var target_body: ?*const types.Body = null;
    var i: usize = 0;
    while (i < simulator.state.bodies.len) : (i += 1) {
        const b = &simulator.state.bodies[i];
        if (std.mem.eql(u8, b.id, target_id)) {
            target_body = b;
            break;
        }
    }

    var player_ship: ?*types.Ship = null;
    i = 0;
    while (i < simulator.state.ships.len) : (i += 1) {
        const s = &simulator.state.ships[i];
        if (std.mem.eql(u8, s.body.id, player_id)) {
            player_ship = s;
            break;
        }
    }

    if (target_body) |t| {
        if (player_ship) |ship| {
            var dir = math.Vec3.sub(ship.body.position, t.position);
            const len = dir.len();
            if (len < 1e-9) {
                dir = types.Vec3{ .x = 1.0, .y = 0.0, .z = 0.0 };
            } else {
                dir = dir.scale(1.0 / len);
            }
            const baseR = if (t.terrain) |terrain| terrain.radius else t.radius;
            const sea = if (t.terrain) |terrain| terrain.seaLevel else 0.0;
            const maxH = if (t.terrain) |terrain| (terrain.maxHeight * baseR) else 0.0;
            const surfaceR = baseR + @max(sea, maxH);
            const MIN_ALT: f64 = 2.1; // km
            const new_pos = math.Vec3.add(t.position, dir.scale(surfaceR + MIN_ALT));
            ship.body.position = new_pos;
            ship.body.velocity = t.velocity;
            ship.angularVelocity = types.Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 };
        }
    }
}

pub export fn auto_land(simulator: *sim.Simulator, target_id_ptr: u32, target_id_len: u32) callconv(.c) void {
    const id_ptr: [*]const u8 = @ptrFromInt(target_id_ptr);
    const target_id = id_ptr[0..target_id_len];
    const player_id = "player-ship";
    // Find target body index
    var target_index: ?usize = null;
    var i: usize = 0;
    while (i < simulator.state.bodies.len) : (i += 1) {
        const b = simulator.state.bodies[i];
        if (std.mem.eql(u8, b.id, target_id)) {
            target_index = i;
            break;
        }
    }
    if (target_index == null) return;

    // Ensure the player ship exists
    var player_exists = false;
    i = 0;
    while (i < simulator.state.ships.len) : (i += 1) {
        const s = simulator.state.ships[i];
        if (std.mem.eql(u8, s.body.id, player_id)) {
            player_exists = true;
            break;
        }
    }
    if (!player_exists) return;

    // Arm suicide-burn autoland controller
    simulator.auto_land_active = true;
    simulator.auto_land_target_index = target_index.?;
}

pub export fn free_result_buffer(ptr_addr: u32, len: u32) callconv(.c) void {
    const allocator = wasmAllocator();
    const base_ptr: [*]u8 = @ptrFromInt(@as(usize, ptr_addr));
    const slice: []u8 = base_ptr[0..len];
    allocator.free(slice);
}

pub export fn destroy_simulator(simulator: *sim.Simulator) callconv(.c) void {
    simulator.destroy();
}

// --- Coarse Grid Bridge ---
var active_grid: ?*planet.Grid = null;

pub export fn create_grid(size: u32) callconv(.c) void {
    const alloc = wasmAllocator();
    if (active_grid) |g| {
        g.deinit();
        alloc.destroy(g);
        active_grid = null;
    }
    const g_ptr = alloc.create(planet.Grid) catch return;
    g_ptr.* = planet.Grid.init(alloc, @intCast(size)) catch {
        alloc.destroy(g_ptr);
        return;
    };
    const g = g_ptr;
    var coord_map = g.populateIndices() catch {
        g.deinit();
        alloc.destroy(g);
        return;
    };
    defer coord_map.deinit();
    g.populateNeighbors(&coord_map) catch {
        g.deinit();
        alloc.destroy(g);
        return;
    };
    g.generateMesh() catch {
        g.deinit();
        alloc.destroy(g);
        return;
    };
    active_grid = g;
}

fn asSliceU8(ptr: [*]const u8, len_bytes: usize) SliceU8 {
    return .{ .ptr = @intCast(@intFromPtr(ptr)), .len = @intCast(len_bytes) };
}

pub export fn get_grid_vertex_buffer() callconv(.c) SliceU8 {
    if (active_grid == null) return .{ .ptr = 0, .len = 0 };
    const g = active_grid.?;
    if (g.vertices == null) return .{ .ptr = 0, .len = 0 };
    const verts = g.vertices.?;
    const bytes: usize = verts.len * @sizeOf(f32);
    const base: [*]const u8 = @ptrCast(verts.ptr);
    return asSliceU8(base, bytes);
}

pub export fn get_grid_vertex_buffer_out(out_ptr: u32) callconv(.c) void {
    const r = get_grid_vertex_buffer();
    const out: [*]u32 = @ptrFromInt(out_ptr);
    out[0] = r.ptr;
    out[1] = r.len;
}

pub export fn get_grid_elevation_buffer() callconv(.c) SliceU8 {
    if (active_grid == null) return .{ .ptr = 0, .len = 0 };
    const g = active_grid.?;
    if (g.elevations == null) return .{ .ptr = 0, .len = 0 };
    const arr = g.elevations.?;
    const bytes: usize = arr.len * @sizeOf(f32);
    const base: [*]const u8 = @ptrCast(arr.ptr);
    return asSliceU8(base, bytes);
}

pub export fn get_grid_elevation_buffer_out(out_ptr: u32) callconv(.c) void {
    const r = get_grid_elevation_buffer();
    const out: [*]u32 = @ptrFromInt(out_ptr);
    out[0] = r.ptr;
    out[1] = r.len;
}

pub export fn get_grid_index_buffer() callconv(.c) SliceU8 {
    if (active_grid == null) return .{ .ptr = 0, .len = 0 };
    const g = active_grid.?;
    if (g.indices == null) return .{ .ptr = 0, .len = 0 };
    const arr = g.indices.?;
    const bytes: usize = arr.len * @sizeOf(u32);
    const base: [*]const u8 = @ptrCast(arr.ptr);
    return asSliceU8(base, bytes);
}

pub export fn get_grid_index_buffer_out(out_ptr: u32) callconv(.c) void {
    const r = get_grid_index_buffer();
    const out: [*]u32 = @ptrFromInt(out_ptr);
    out[0] = r.ptr;
    out[1] = r.len;
}

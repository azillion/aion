const std = @import("std");
const json = std.json;
const types = @import("types.zig");

// Import JS logging function for error reporting
extern "env" fn log_error(ptr: [*]const u8, len: usize) void;

const G: f64 = 6.67430e-20; // km^3 / (kg * s^2)

// C-compatible slice return for WASM FFI
pub const SliceU8 = extern struct { ptr: u32, len: u32 };

pub const Simulator = struct {
    allocator: std.mem.Allocator,
    state: types.SystemState,

    pub fn create(allocator: std.mem.Allocator, initial_state: types.SystemState) !*Simulator {
        const sim = try allocator.create(Simulator);
        sim.* = .{ .allocator = allocator, .state = initial_state };
        return sim;
    }

    pub fn destroy(self: *Simulator) void {
        // free string fields for bodies then free slice
        for (self.state.bodies) |b| {
            self.allocator.free(b.id);
            self.allocator.free(b.name);
        }
        self.allocator.free(self.state.bodies);

        // free string fields for ships then free slice
        for (self.state.ships) |s| {
            self.allocator.free(s.id);
            self.allocator.free(s.name);
        }
        self.allocator.free(self.state.ships);

        self.allocator.destroy(self);
    }
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

pub export fn get_input_buffer() callconv(.c) u32 {
    return @intCast(@intFromPtr(&input_buffer));
}

pub export fn create_simulator(initial_state_ptr: u32, initial_state_len: u32) callconv(.c) ?*Simulator {
    const allocator = wasmAllocator();
    const ptr: [*]const u8 = @ptrFromInt(initial_state_ptr);
    const input_slice = ptr[0..initial_state_len];
    var parsed = json.parseFromSlice(types.SystemState, allocator, input_slice, .{}) catch |err| {
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
    const sim = Simulator.create(allocator, owned_state) catch {
        parsed.deinit();
        return null;
    };
    // Intentionally do not parsed.deinit(); simulator owns the memory
    return sim;
}

fn quatNormalize(q: *types.Quat) void {
    const n = @sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (n > 0) {
        q[0] /= n;
        q[1] /= n;
        q[2] /= n;
        q[3] /= n;
    } else {
        q.* = .{ 0, 0, 0, 1 };
    }
}

fn quatMultiply(a: types.Quat, b: types.Quat) types.Quat {
    return types.Quat{
        a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
        a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
        a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
        a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2],
    };
}

fn quatFromAxisAngle(axis: types.Vec3, angle: f64) types.Quat {
    const half = angle * 0.5;
    const s = @sin(half);
    const c = @cos(half);
    return types.Quat{ axis[0] * s, axis[1] * s, axis[2] * s, c };
}

fn rotateVecByQuat(v: types.Vec3, q: types.Quat) types.Vec3 {
    // v' = v + 2*cross(q.xyz, cross(q.xyz, v) + q.w*v)
    const ux = q[0];
    const uy = q[1];
    const uz = q[2];
    const s = q[3];
    const vx = v[0];
    const vy = v[1];
    const vz = v[2];
    const cx = uy * vz - uz * vy;
    const cy = uz * vx - ux * vz;
    const cz = ux * vy - uy * vx;
    const c2x = uy * cz - uz * cy + s * vx;
    const c2y = uz * cx - ux * cz + s * vy;
    const c2z = ux * cy - uy * cx + s * vz;
    return .{ vx + 2.0 * (uy * c2z - uz * c2y), vy + 2.0 * (uz * c2x - ux * c2z), vz + 2.0 * (ux * c2y - uy * c2x) };
}

fn vec3Len(v: [3]f64) f64 {
    return @sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
fn vec3Scale(v: [3]f64, s: f64) [3]f64 {
    return .{ v[0] * s, v[1] * s, v[2] * s };
}
fn vec3Add(a: [3]f64, b: [3]f64) [3]f64 {
    return .{ a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

pub export fn tick_simulator(sim: *Simulator, dt: f64) callconv(.c) void {
    const bodies = sim.state.bodies;
    const ships = sim.state.ships;
    const num_bodies = bodies.len;
    const num_ships = ships.len;
    const n = num_bodies + num_ships;
    if (n == 0) return;

    var accelerations = sim.allocator.alloc(types.Vec3, n) catch return;
    defer sim.allocator.free(accelerations);
    for (accelerations) |*acc| acc.* = .{ 0.0, 0.0, 0.0 };

    // Helpers to read/write by combined index
    const All = struct {
        bodies: []const types.Body,
        ships: []const types.Ship,
        fn pos(self: @This(), idx: usize) types.Vec3 {
            return if (idx < self.bodies.len) self.bodies[idx].position else self.ships[idx - self.bodies.len].position;
        }
        fn mass(self: @This(), idx: usize) f64 {
            return if (idx < self.bodies.len) self.bodies[idx].mass else self.ships[idx - self.bodies.len].mass;
        }
    };
    const all = All{ .bodies = bodies, .ships = ships };

    // Pairwise gravity across bodies + ships
    var i: usize = 0;
    while (i < n) : (i += 1) {
        var j: usize = i + 1;
        while (j < n) : (j += 1) {
            const a_pos = all.pos(i);
            const b_pos = all.pos(j);
            const dx = b_pos[0] - a_pos[0];
            const dy = b_pos[1] - a_pos[1];
            const dz = b_pos[2] - a_pos[2];
            const dist_sq = dx * dx + dy * dy + dz * dz;
            const safe = if (dist_sq < 1e-9) 1e-9 else dist_sq;
            const dist = @sqrt(safe);
            const force_mag = (G * all.mass(i) * all.mass(j)) / safe;
            const fx = force_mag * (dx / dist);
            const fy = force_mag * (dy / dist);
            const fz = force_mag * (dz / dist);
            accelerations[i][0] += fx / all.mass(i);
            accelerations[i][1] += fy / all.mass(i);
            accelerations[i][2] += fz / all.mass(i);
            accelerations[j][0] -= fx / all.mass(j);
            accelerations[j][1] -= fy / all.mass(j);
            accelerations[j][2] -= fz / all.mass(j);
        }
    }

    // Ship controls (only ships)
    var player_ship_index: ?usize = null;
    {
        const target_id = "player-ship";
        var s: usize = 0;
        while (s < num_ships) : (s += 1) {
            if (std.mem.eql(u8, ships[s].id, target_id)) {
                player_ship_index = s;
                break;
            }
        }
    }
    if (player_ship_index) |si| {
        const ship = &sim.state.ships[si];
        const acc_index = num_bodies + si;
        const mask = input_buffer.keys_mask;
        const isSet = struct {
            fn f(m: u32, bit: u5) bool {
                return (m & (@as(u32, 1) << bit)) != 0;
            }
        };
        const PREC = isSet.f(mask, 14);
        const precision_multiplier: f64 = if (PREC) 0.1 else 1.0;
        const THRUST_FORCE: f64 = 5e5 * precision_multiplier;
        const MANEUVER_FORCE: f64 = THRUST_FORCE / 4.0;
        const ANGULAR_ACCELERATION: f64 = 0.5 * precision_multiplier;

        const base_forward: types.Vec3 = .{ 0, 0, -1 };
        const base_right: types.Vec3 = .{ 1, 0, 0 };
        const base_up: types.Vec3 = .{ 0, 1, 0 };
        const q: types.Quat = ship.orientation;
        const ship_forward = rotateVecByQuat(base_forward, q);
        const ship_right = rotateVecByQuat(base_right, q);
        const ship_up = rotateVecByQuat(base_up, q);

        var total_force: [3]f64 = .{ 0, 0, 0 };
        if (isSet.f(mask, 0)) total_force = vec3Add(total_force, vec3Scale(ship_forward, THRUST_FORCE));
        if (isSet.f(mask, 1)) total_force = vec3Add(total_force, vec3Scale(ship_forward, -THRUST_FORCE / 2.0));
        if (isSet.f(mask, 3)) total_force = vec3Add(total_force, vec3Scale(ship_right, MANEUVER_FORCE));
        if (isSet.f(mask, 2)) total_force = vec3Add(total_force, vec3Scale(ship_right, -MANEUVER_FORCE));
        if (isSet.f(mask, 4)) total_force = vec3Add(total_force, vec3Scale(ship_up, MANEUVER_FORCE));
        if (isSet.f(mask, 5)) total_force = vec3Add(total_force, vec3Scale(ship_up, -MANEUVER_FORCE));

        var ang_vel: [3]f64 = ship.angularVelocity;
        if (isSet.f(mask, 6)) {
            const DAMP = 0.95;
            ang_vel = vec3Scale(ang_vel, DAMP);
            if (vec3Len(ang_vel) < 0.001) ang_vel = .{ 0, 0, 0 };
        } else {
            var torque: [3]f64 = .{ 0, 0, 0 };
            if (isSet.f(mask, 7)) torque[0] -= ANGULAR_ACCELERATION;
            if (isSet.f(mask, 8)) torque[0] += ANGULAR_ACCELERATION;
            if (isSet.f(mask, 9)) torque[1] += ANGULAR_ACCELERATION;
            if (isSet.f(mask, 10)) torque[1] -= ANGULAR_ACCELERATION;
            if (isSet.f(mask, 11)) torque[2] += ANGULAR_ACCELERATION;
            if (isSet.f(mask, 12)) torque[2] -= ANGULAR_ACCELERATION;
            ang_vel = vec3Add(ang_vel, vec3Scale(torque, dt));
        }
        ship.angularVelocity = ang_vel;

        const omega = ang_vel;
        const angle = vec3Len(omega) * dt;
        if (angle > 0) {
            const axis = vec3Scale(omega, 1.0 / @max(vec3Len(omega), 1e-9));
            const dq: types.Quat = quatFromAxisAngle(axis, angle);
            var qcur: types.Quat = ship.orientation;
            qcur = quatMultiply(qcur, dq);
            ship.orientation = qcur;
            quatNormalize(&ship.orientation);
        }

        // Braking (X: bit 13)
        if (isSet.f(mask, 13)) {
            const v = ship.velocity;
            const v_forward = v[0] * ship_forward[0] + v[1] * ship_forward[1] + v[2] * ship_forward[2];
            if (@abs(v_forward) > 0.1) {
                const brake_dir = vec3Scale(ship_forward, if (v_forward > 0) -1.0 else 1.0);
                total_force = vec3Add(total_force, vec3Scale(brake_dir, THRUST_FORCE));
            }
        }

        const acc = vec3Scale(total_force, 1.0 / @max(ship.mass, 1e-9));
        accelerations[acc_index][0] += acc[0];
        accelerations[acc_index][1] += acc[1];
        accelerations[acc_index][2] += acc[2];
    }

    // Integrate both bodies and ships
    i = 0;
    while (i < num_bodies) : (i += 1) {
        var b = &sim.state.bodies[i];
        b.velocity[0] += accelerations[i][0] * dt;
        b.velocity[1] += accelerations[i][1] * dt;
        b.velocity[2] += accelerations[i][2] * dt;
        b.position[0] += b.velocity[0] * dt;
        b.position[1] += b.velocity[1] * dt;
        b.position[2] += b.velocity[2] * dt;
    }
    i = 0;
    while (i < num_ships) : (i += 1) {
        const ai = num_bodies + i;
        var s = &sim.state.ships[i];
        s.velocity[0] += accelerations[ai][0] * dt;
        s.velocity[1] += accelerations[ai][1] * dt;
        s.velocity[2] += accelerations[ai][2] * dt;
        s.position[0] += s.velocity[0] * dt;
        s.position[1] += s.velocity[1] * dt;
        s.position[2] += s.velocity[2] * dt;
    }
}

pub export fn query_simulator_state(sim: *Simulator) SliceU8 {
    const alloc = wasmAllocator();
    // Use ArrayListUnmanaged; provide allocator at usage sites
    var a = std.ArrayListUnmanaged(u8){};
    defer a.deinit(alloc);

    var old = a.writer(alloc);
    var ad = old.adaptToNewApi(&.{});
    const w: *std.Io.Writer = &ad.new_interface;

    const fmt = std.json.fmt(sim.state, .{});
    fmt.format(w) catch return .{ .ptr = 0, .len = 0 };

    const buf = a.toOwnedSlice(alloc) catch return .{ .ptr = 0, .len = 0 };
    return .{ .ptr = @intCast(@intFromPtr(buf.ptr)), .len = @intCast(buf.len) };
}

// Out-parameter variant to avoid ABI differences across engines.
pub export fn query_simulator_state_out(sim: *Simulator, out_ptr: u32) callconv(.c) void {
    const alloc = wasmAllocator();

    var a = std.ArrayListUnmanaged(u8){};
    defer a.deinit(alloc);

    var old = a.writer(alloc);
    var ad = old.adaptToNewApi(&.{});
    const w: *std.Io.Writer = &ad.new_interface;

    const fmt = std.json.fmt(sim.state, .{});
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

pub export fn free_result_buffer(ptr_addr: u32, len: u32) callconv(.c) void {
    const allocator = wasmAllocator();
    const base_ptr: [*]u8 = @ptrFromInt(@as(usize, ptr_addr));
    const slice: []u8 = base_ptr[0..len];
    allocator.free(slice);
}

pub export fn destroy_simulator(sim: *Simulator) callconv(.c) void {
    sim.destroy();
}

const std = @import("std");
const json = std.json;
const types = @import("types.zig");

const G: f64 = 6.67430e-20; // km^3 / (kg * s^2)

pub const Simulator = struct {
    allocator: std.mem.Allocator,
    state: types.SystemState,

    pub fn create(allocator: std.mem.Allocator, initial_state: types.SystemState) !*Simulator {
        const sim = try allocator.create(Simulator);
        sim.* = .{ .allocator = allocator, .state = initial_state };
        return sim;
    }

    pub fn destroy(self: *Simulator) void {
        self.state.bodies.deinit(self.allocator);
        self.allocator.destroy(self);
    }
};

fn wasmAllocator() std.mem.Allocator {
    const builtin = @import("builtin");
    return if (builtin.target.cpu.arch == .wasm32) std.heap.wasm_allocator else std.heap.page_allocator;
}

export fn create_simulator(initial_state_ptr: [*]const u8, initial_state_len: usize) callconv(.c) ?*Simulator {
    const allocator = wasmAllocator();
    const input_slice = initial_state_ptr[0..initial_state_len];
    var parsed = json.parseFromSlice(types.SystemState, allocator, input_slice, .{}) catch return null;
    defer parsed.deinit();

    var owned_state = types.SystemState{
        .timestamp = parsed.value.timestamp,
        .bodies = parsed.value.bodies.clone(allocator) catch return null,
    };

    const sim = Simulator.create(allocator, owned_state) catch {
        owned_state.bodies.deinit(allocator);
        return null;
    };
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

export fn tick_simulator(sim: *Simulator, dt: f64, input_ptr: *const types.InputState) callconv(.c) void {
    var bodies = sim.state.bodies.items;
    const n = bodies.len;
    if (n == 0) return;

    var accelerations = sim.allocator.alloc(types.Vec3, n) catch return;
    defer sim.allocator.free(accelerations);

    // zero out accelerations
    for (accelerations) |*a| {
        a.* = .{ 0.0, 0.0, 0.0 };
    }

    // pairwise gravity
    var i: usize = 0;
    while (i < n) : (i += 1) {
        var j: usize = i + 1;
        while (j < n) : (j += 1) {
            const bodyA = &bodies[i];
            const bodyB = &bodies[j];

            const dx = bodyB.position[0] - bodyA.position[0];
            const dy = bodyB.position[1] - bodyA.position[1];
            const dz = bodyB.position[2] - bodyA.position[2];
            const dist_sq = dx * dx + dy * dy + dz * dz;
            const safe = if (dist_sq < 1e-9) 1e-9 else dist_sq;
            const dist = @sqrt(safe);
            const force_mag = (G * bodyA.mass * bodyB.mass) / safe;
            const fx = force_mag * (dx / dist);
            const fy = force_mag * (dy / dist);
            const fz = force_mag * (dz / dist);

            accelerations[i][0] += fx / bodyA.mass;
            accelerations[i][1] += fy / bodyA.mass;
            accelerations[i][2] += fz / bodyA.mass;

            accelerations[j][0] -= fx / bodyB.mass;
            accelerations[j][1] -= fy / bodyB.mass;
            accelerations[j][2] -= fz / bodyB.mass;
        }
    }

    // Ship controls (player-ship)
    var player_index: ?usize = null;
    {
        const target_id = "player-ship";
        var k: usize = 0;
        while (k < n) : (k += 1) {
            if (std.mem.eql(u8, bodies[k].id, target_id)) {
                player_index = k;
                break;
            }
        }
    }
    if (player_index) |pi| {
        const ship = &bodies[pi];
        const mask = input_ptr.keys_mask;
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

        // Basis in ship local space
        const base_forward: types.Vec3 = .{ 0, 0, -1 };
        const base_right: types.Vec3 = .{ 1, 0, 0 };
        const base_up: types.Vec3 = .{ 0, 1, 0 };
        const q: types.Quat = ship.orientation;
        const ship_forward = rotateVecByQuat(base_forward, q);
        const ship_right = rotateVecByQuat(base_right, q);
        const ship_up = rotateVecByQuat(base_up, q);

        var total_force: [3]f64 = .{ 0, 0, 0 };
        if (isSet.f(mask, 0)) {
            total_force = vec3Add(total_force, vec3Scale(ship_forward, THRUST_FORCE));
        }
        if (isSet.f(mask, 1)) {
            total_force = vec3Add(total_force, vec3Scale(ship_forward, -THRUST_FORCE / 2.0));
        }
        if (isSet.f(mask, 3)) {
            total_force = vec3Add(total_force, vec3Scale(ship_right, MANEUVER_FORCE));
        }
        if (isSet.f(mask, 2)) {
            total_force = vec3Add(total_force, vec3Scale(ship_right, -MANEUVER_FORCE));
        }
        if (isSet.f(mask, 4)) {
            total_force = vec3Add(total_force, vec3Scale(ship_up, MANEUVER_FORCE));
        }
        if (isSet.f(mask, 5)) {
            total_force = vec3Add(total_force, vec3Scale(ship_up, -MANEUVER_FORCE));
        }

        // Angular controls
        var ang_vel: [3]f64 = ship.angularVelocity;
        if (isSet.f(mask, 6)) { // Backspace: dampen
            const DAMP = 0.95;
            ang_vel = vec3Scale(ang_vel, DAMP);
            if (vec3Len(ang_vel) < 0.001) ang_vel = .{ 0, 0, 0 };
        } else {
            var torque: [3]f64 = .{ 0, 0, 0 };
            // Pitch R/F (indices 7,8)
            if (isSet.f(mask, 7)) {
                torque[0] -= ANGULAR_ACCELERATION;
            }
            if (isSet.f(mask, 8)) {
                torque[0] += ANGULAR_ACCELERATION;
            }
            // Yaw Z/C (9,10)
            if (isSet.f(mask, 9)) {
                torque[1] += ANGULAR_ACCELERATION;
            }
            if (isSet.f(mask, 10)) {
                torque[1] -= ANGULAR_ACCELERATION;
            }
            // Roll Q/E (11,12)
            if (isSet.f(mask, 11)) {
                torque[2] += ANGULAR_ACCELERATION;
            }
            if (isSet.f(mask, 12)) {
                torque[2] -= ANGULAR_ACCELERATION;
            }
            ang_vel = vec3Add(ang_vel, vec3Scale(torque, dt));
        }
        ship.angularVelocity = ang_vel;

        // Orientation integrate from angular velocity
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

        // Add thrust acceleration to accelerations
        const acc = vec3Scale(total_force, 1.0 / @max(ship.mass, 1e-9));
        accelerations[pi][0] += acc[0];
        accelerations[pi][1] += acc[1];
        accelerations[pi][2] += acc[2];
    }

    // integrate (Euler)
    i = 0;
    while (i < n) : (i += 1) {
        var body = &bodies[i];
        body.velocity[0] += accelerations[i][0] * dt;
        body.velocity[1] += accelerations[i][1] * dt;
        body.velocity[2] += accelerations[i][2] * dt;

        body.position[0] += body.velocity[0] * dt;
        body.position[1] += body.velocity[1] * dt;
        body.position[2] += body.velocity[2] * dt;
    }
}

export fn query_simulator_state(sim: *Simulator) callconv(.c) u64 {
    var w: std.Io.Writer.Allocating = .init(sim.allocator);
    std.json.Stringify.value(sim.state, .{}, &w.writer) catch return 0;
    const buf = w.toOwnedSlice() catch return 0;
    const addr: u64 = @intFromPtr(buf.ptr);
    const len: u64 = buf.len;
    return (addr << 32) | len;
}

export fn free_result_buffer(ptr_addr: u32, len: u32) callconv(.c) void {
    const allocator = wasmAllocator();
    const base_ptr: [*]u8 = @ptrFromInt(@as(usize, ptr_addr));
    const slice: []u8 = base_ptr[0..len];
    allocator.free(slice);
}

export fn destroy_simulator(sim: *Simulator) callconv(.c) void {
    sim.destroy();
}

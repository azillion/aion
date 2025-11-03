const std = @import("std");
const std_math = @import("std").math;
const types = @import("core");
const math = @import("math");

const G: f64 = 6.67430e-20; // km^3 / (kg * s^2)

pub fn tickSimulator(simulator: *anyopaque, dt_unscaled: f64, input_buffer: *types.InputState) void {
    // Cast simulator pointer to the actual type
    const sim = @as(*@import("sim").Simulator, @ptrCast(@alignCast(simulator)));

    const dt = dt_unscaled * sim.time_scale;
    if (dt == 0.0) return;
    const bodies = sim.state.bodies;
    const ships = sim.state.ships;
    const num_bodies = bodies.len;
    const num_ships = ships.len;
    const n = num_bodies + num_ships;
    if (n == 0) return;

    var accelerations = sim.allocator.alloc(types.Vec3, n) catch return;
    defer sim.allocator.free(accelerations);
    for (accelerations) |*acc| acc.* = math.Vec3.ZERO;

    // Helpers to read/write by combined index
    const All = struct {
        bodies: []const types.Body,
        ships: []const types.Ship,
        fn pos(self: @This(), idx: usize) types.Vec3 {
            return if (idx < self.bodies.len) self.bodies[idx].position else self.ships[idx - self.bodies.len].body.position;
        }
        fn mass(self: @This(), idx: usize) f64 {
            return if (idx < self.bodies.len) self.bodies[idx].mass else self.ships[idx - self.bodies.len].body.mass;
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
            const diff = math.Vec3.sub(b_pos, a_pos);
            const dist_sq = math.Vec3.len2(diff);
            const safe = if (dist_sq < 1e-9) 1e-9 else dist_sq;
            const dist = @sqrt(safe);
            const force_mag = (G * all.mass(i) * all.mass(j)) / safe;
            const force_dir = diff.scale(1.0 / dist);
            const force = force_dir.scale(force_mag);
            accelerations[i] = math.Vec3.add(accelerations[i], force.scale(1.0 / all.mass(i)));
            accelerations[j] = math.Vec3.sub(accelerations[j], force.scale(1.0 / all.mass(j)));
        }
    }

    // Ship controls (only ships)
    var player_ship_index: ?usize = null;
    {
        const target_id = "player-ship";
        var s: usize = 0;
        while (s < num_ships) : (s += 1) {
            if (std.mem.eql(u8, ships[s].body.id, target_id)) {
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

        const base_forward = math.Vec3.UNIT_Z.negate();
        const base_right = math.Vec3.UNIT_X;
        const base_up = math.Vec3.UNIT_Y;
        const q: types.Quat = ship.orientation;
        const ship_forward = math.Quat.rotateVec(base_forward, q);
        const ship_right = math.Quat.rotateVec(base_right, q);
        const ship_up = math.Quat.rotateVec(base_up, q);

        var total_force = math.Vec3.ZERO;
        if (isSet.f(mask, 0)) total_force = math.Vec3.add(total_force, ship_forward.scale(THRUST_FORCE));
        if (isSet.f(mask, 1)) total_force = math.Vec3.add(total_force, ship_forward.scale(-THRUST_FORCE / 2.0));
        if (isSet.f(mask, 3)) total_force = math.Vec3.add(total_force, ship_right.scale(MANEUVER_FORCE));
        if (isSet.f(mask, 2)) total_force = math.Vec3.add(total_force, ship_right.scale(-MANEUVER_FORCE));
        if (isSet.f(mask, 4)) total_force = math.Vec3.add(total_force, ship_up.scale(MANEUVER_FORCE));
        if (isSet.f(mask, 5)) total_force = math.Vec3.add(total_force, ship_up.scale(-MANEUVER_FORCE));

        var ang_vel = ship.angularVelocity;
        // Flight assist (debug): KeyB on bit 15 instantly snaps orientation to face Sun
        if (isSet.f(mask, 15)) {
            const sun_id = "sol";
            var sun: ?*const types.Body = null;
            var bi: usize = 0;
            while (bi < sim.state.bodies.len) : (bi += 1) {
                const b = &sim.state.bodies[bi];
                if (std.mem.eql(u8, b.id, sun_id)) {
                    sun = b;
                    break;
                }
            }
            if (sun) |s| {
                const target_vec = math.Vec3.sub(s.position, ship.body.position);
                ship.orientation = math.Quat.fromTwoVectors(base_forward, target_vec);
                math.Quat.normalize(&ship.orientation);
                ang_vel = types.Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 };
            }
        } else if (isSet.f(mask, 6)) {
            const DAMP = 0.95;
            ang_vel = ang_vel.scale(DAMP);
            if (ang_vel.len() < 0.001) ang_vel = math.Vec3.ZERO;
        } else {
            var torque = math.Vec3.ZERO;
            if (isSet.f(mask, 7)) torque.x -= ANGULAR_ACCELERATION;
            if (isSet.f(mask, 8)) torque.x += ANGULAR_ACCELERATION;
            if (isSet.f(mask, 9)) torque.y += ANGULAR_ACCELERATION;
            if (isSet.f(mask, 10)) torque.y -= ANGULAR_ACCELERATION;
            if (isSet.f(mask, 11)) torque.z += ANGULAR_ACCELERATION;
            if (isSet.f(mask, 12)) torque.z -= ANGULAR_ACCELERATION;
            ang_vel = math.Vec3.add(ang_vel, torque.scale(dt));
        }
        ship.angularVelocity = ang_vel;

        const omega = ang_vel;
        const angle = omega.len() * dt;
        if (angle > 0) {
            const axis = omega.scale(1.0 / @max(omega.len(), 1e-9));
            const dq = math.Quat.fromAxisAngle(axis, angle);
            ship.orientation = math.Quat.multiply(ship.orientation, dq);
            math.Quat.normalize(&ship.orientation);
        }

        // Braking (X: bit 13)
        if (isSet.f(mask, 13)) {
            const v = ship.body.velocity;
            const v_forward = math.Vec3.dot(v, ship_forward);
            if (@abs(v_forward) > 0.1) {
                const brake_dir = ship_forward.scale(if (v_forward > 0) -1.0 else 1.0);
                total_force = math.Vec3.add(total_force, brake_dir.scale(THRUST_FORCE));
            }
        }

        const acc = total_force.scale(1.0 / @max(ship.body.mass, 1e-9));
        accelerations[acc_index] = math.Vec3.add(accelerations[acc_index], acc);

        // --- Advanced autoland: gravity-aware bang-bang with tangential cancellation ---
        if (sim.auto_land_active) {
            const target = sim.state.bodies[sim.auto_land_target_index];
            var diff_in = math.Vec3.sub(target.position, ship.body.position); // inward toward target
            const dist = diff_in.len();
            if (dist > 0.0) {
                const n_in = diff_in.scale(1.0 / dist);
                const n_out = n_in.scale(-1.0);
                const baseR = if (target.terrain) |t| t.radius else target.radius;
                const sea = if (target.terrain) |t| t.seaLevel else 0.0;
                const maxH = if (target.terrain) |t| (t.maxHeight * baseR) else 0.0;
                const surfaceR = baseR + @max(sea, maxH);
                const minAlt = sim.auto_land_min_alt_km;
                const minDist = surfaceR + minAlt;

                const v_rel = math.Vec3.sub(ship.body.velocity, target.velocity);
                const v_in = v_rel.x * n_in.x + v_rel.y * n_in.y + v_rel.z * n_in.z; // inward +
                const v_in_vec = n_in.scale(v_in);
                const v_tan = math.Vec3.sub(v_rel, v_in_vec);
                const v_tan_len = v_tan.len();
                const r = dist - minDist;

                // Thrust capability
                const AUTOPILOT_FORCE: f64 = 5e5; // matches THRUST_FORCE scale
                const a_max = AUTOPILOT_FORCE / @max(ship.body.mass, 1e-9);

                // Estimate gravity toward target for switching rule
                const mu = 6.67430e-20 * target.mass; // km^3/s^2
                const g_in = mu / @max(dist * dist, 1e-9); // inward (+)

                // Allocate some acceleration to cancel tangential velocity
                // Aim to kill tangential over ~2 seconds, capped to 60% of a_max
                var a_tan_mag = v_tan_len / 2.0;
                if (a_tan_mag > a_max * 0.6) a_tan_mag = a_max * 0.6;
                const a_rad_budget = @max(a_max - a_tan_mag, 0.0);

                var a_cmd = math.Vec3.ZERO;

                if (r <= 0.0) {
                    // Inside target altitude: clamp and finish
                    ship.body.position = math.Vec3.add(target.position, n_out.scale(minDist));
                    ship.body.velocity = target.velocity;
                    ship.angularVelocity = types.Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 };
                    sim.auto_land_active = false;
                } else {
                    // Gravity-aware brake distance
                    const denom = 2.0 * (a_rad_budget + g_in + 1e-9);
                    const brakeDist = if (v_in > 0.0) ((v_in * v_in) / denom) else 0.0;

                    // Radial command: inward accelerate until switching surface, then outward brake
                    if ((v_in > 0.0) and (r <= brakeDist + 0.5)) {
                        // Brake outward using full radial budget
                        a_cmd = math.Vec3.add(a_cmd, n_out.scale(a_rad_budget));
                    } else if ((v_in <= 0.0) or (r > brakeDist)) {
                        // Accelerate inward to build approach speed
                        a_cmd = math.Vec3.add(a_cmd, n_in.scale(a_rad_budget));
                    }

                    // Tangential cancellation using remaining budget slice
                    if (v_tan_len > 0.01 and a_tan_mag > 0.0) {
                        const tan_dir = v_tan.scale(1.0 / @max(v_tan_len, 1e-9));
                        a_cmd = math.Vec3.add(a_cmd, tan_dir.scale(-a_tan_mag));
                    }

                    // Cap total magnitude to a_max
                    const a_len = a_cmd.len();
                    if (a_len > a_max) {
                        a_cmd = a_cmd.scale(a_max / a_len);
                    }

                    // Terminal settle: close and slow
                    const speed_rel = v_rel.len();
                    if ((r < 0.5) and (speed_rel < 0.05)) {
                        ship.body.position = math.Vec3.add(target.position, n_out.scale(minDist));
                        ship.body.velocity = target.velocity;
                        ship.angularVelocity = types.Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 };
                        sim.auto_land_active = false;
                        a_cmd = math.Vec3.ZERO;
                    }

                    // Apply commanded acceleration
                    accelerations[acc_index] = math.Vec3.add(accelerations[acc_index], a_cmd);

                    // Point ship along thrust direction when burning
                    if (a_cmd.len() > 0.0) {
                        const desired_dir = a_cmd.scale(1.0 / @max(a_cmd.len(), 1e-9));
                        ship.orientation = math.Quat.fromTwoVectors(base_forward, desired_dir);
                        math.Quat.normalize(&ship.orientation);
                        ship.angularVelocity = types.Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 };
                    }
                }
            }
        }
    }

    // Integrate both bodies and ships
    i = 0;
    while (i < num_bodies) : (i += 1) {
        var b = &sim.state.bodies[i];
        b.velocity = math.Vec3.add(b.velocity, accelerations[i].scale(dt));
        b.position = math.Vec3.add(b.position, b.velocity.scale(dt));
    }
    i = 0;
    while (i < num_ships) : (i += 1) {
        const ai = num_bodies + i;
        var s = &sim.state.ships[i];
        s.body.velocity = math.Vec3.add(s.body.velocity, accelerations[ai].scale(dt));
        s.body.position = math.Vec3.add(s.body.position, s.body.velocity.scale(dt));

        // --- Simple ground clamp at 2 km AGL relative to nearest massive body ---
        var nearest_index: ?usize = null;
        var nearest_dist: f64 = 1.0e300;
        var j: usize = 0;
        while (j < num_bodies) : (j += 1) {
            const b = sim.state.bodies[j];
            if (b.mass < 1e22) continue; // only clamp to sufficiently massive bodies
            const diff = math.Vec3.sub(s.body.position, b.position);
            const d2 = math.Vec3.len2(diff);
            if (d2 < nearest_dist) {
                nearest_dist = d2;
                nearest_index = j;
            }
        }
        if (nearest_index) |ni| {
            const b = sim.state.bodies[ni];
            const diff = math.Vec3.sub(s.body.position, b.position);
            const dist = math.Vec3.len(diff);
            if (dist > 0) {
                // Compute physical surface radius (km)
                const baseR = if (b.terrain) |t| t.radius else b.radius;
                const sea = if (b.terrain) |t| t.seaLevel else 0.0;
                const maxH = if (b.terrain) |t| (t.maxHeight * baseR) else 0.0;
                const surfaceR = baseR + @max(sea, maxH);
                const minAltitude = 2.0; // km AGL clamp
                const minDist = surfaceR + minAltitude;
                if (dist < minDist) {
                    // Move ship out to minDist along the outward normal
                    const normal = diff.scale(1.0 / dist);
                    s.body.position = math.Vec3.add(b.position, normal.scale(minDist));
                    // Hard hold: match planet velocity at clamp altitude to avoid energy injection
                    s.body.velocity = b.velocity;
                }
            }
        }
    }
}

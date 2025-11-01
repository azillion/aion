const std = @import("std");
const vec3 = @import("vec3.zig");
const Vec3 = vec3.Vec3;

pub const Quat = struct {
    const Self = @This();

    x: f64,
    y: f64,
    z: f64,
    w: f64,

    pub const unit = Quat{ .x = 0, .y = 0, .z = 0, .w = 1 };
    pub const zero = Quat{ .x = 0, .y = 0, .z = 0, .w = 0 };

    pub fn normalize(q: *Quat) void {
        const n = @sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
        if (n > 0) {
            q.x /= n;
            q.y /= n;
            q.z /= n;
            q.w /= n;
        } else {
            q.* = Quat{ .x = 0, .y = 0, .z = 0, .w = 1 };
        }
    }

    pub fn multiply(a: Quat, b: Quat) Quat {
        return Quat{
            .x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            .y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            .z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
            .w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        };
    }

    pub fn fromAxisAngle(axis: Vec3, angle: f64) Quat {
        const half = angle * 0.5;
        const s = @sin(half);
        const c = @cos(half);
        return Quat{ .x = axis.x * s, .y = axis.y * s, .z = axis.z * s, .w = c };
    }

    pub fn fromDirection(normal: Vec3, up: ?Vec3) Quat {
        const u = up orelse Vec3{ .x = 0.0, .y = 0.0, .z = 1.0 };
        const n = normal.normalize();
        const a = Vec3.cross(u, n);
        const d = Vec3.dot(u, n);
        return Quat{ .x = a.x, .y = a.y, .z = a.z, .w = d + 1.0 };
    }

    pub fn fromTwoVectors(a: Vec3, b: Vec3) Quat {
        const len_a = a.len();
        const len_b = b.len();
        const norm_a = if (len_a > 1e-9) a.scale(1.0 / len_a) else Vec3{ .x = 0.0, .y = 0.0, .z = -1.0 };
        const norm_b = if (len_b > 1e-9) b.scale(1.0 / len_b) else Vec3{ .x = 0.0, .y = 0.0, .z = -1.0 };

        const dot_prod = norm_a.x * norm_b.x + norm_a.y * norm_b.y + norm_a.z * norm_b.z;
        if (dot_prod < -0.999999) {
            var axis = Vec3{ .x = 1.0, .y = 0.0, .z = 0.0 };
            if (std.math.fabs(norm_a.x) > 0.9) axis = Vec3{ .x = 0.0, .y = 1.0, .z = 0.0 };
            const ortho = vec3.Vec3.cross(axis, norm_a);
            return fromAxisAngle(ortho, std.math.pi);
        }

        const cross = vec3.Vec3.cross(norm_a, norm_b);
        const s = @sqrt((1.0 + dot_prod) * 2.0);
        const invs = 1.0 / @max(s, 1e-9);
        return Quat{ .x = cross.x * invs, .y = cross.y * invs, .z = cross.z * invs, .w = s * 0.5 };
    }

    pub fn add(a: Quat, b: Quat) Quat {
        return Quat{
            .x = a.x + b.x,
            .y = a.y + b.y,
            .z = a.z + b.z,
            .w = a.w + b.w,
        };
    }

    pub fn sub(a: Quat, b: Quat) Quat {
        return Quat{
            .x = a.x - b.x,
            .y = a.y - b.y,
            .z = a.z - b.z,
            .w = a.w - b.w,
        };
    }

    pub fn scale(a: Quat, s: f64) Quat {
        return Quat{
            .x = a.x * s,
            .y = a.y * s,
            .z = a.z * s,
            .w = a.w * s,
        };
    }

    pub fn dot(a: Quat, b: Quat) f64 {
        return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    }

    pub fn len(a: Quat) f64 {
        return @sqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w);
    }

    pub fn len2(a: Quat) f64 {
        return a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w;
    }

    pub fn conjugate(a: Quat) Quat {
        return Quat{
            .x = -a.x,
            .y = -a.y,
            .z = -a.z,
            .w = a.w,
        };
    }

    pub fn inverse(a: Quat) Quat {
        const conj = conjugate(a);
        const len_sq = len2(a);
        if (len_sq > 0) {
            return scale(conj, 1.0 / len_sq);
        } else {
            return zero;
        }
    }

    pub fn slerp(a: Quat, b: Quat, s: f64) Quat {
        var dot_val = dot(a, b);
        var a_mod = a;

        if (dot_val < 0) {
            a_mod = scale(a, -1.0);
            dot_val = -dot_val;
        }

        if (dot_val > 0.9995) {
            const result = add(a_mod, scale(sub(b, a_mod), s));
            var normalized = result;
            normalize(&normalized);
            return normalized;
        }

        const clamped_dot = @max(-1.0, @min(1.0, dot_val));
        const theta = std.math.acos(clamped_dot) * s;

        const c = normalize(sub(b, scale(a_mod, dot_val)));
        return add(scale(a_mod, @cos(theta)), scale(c, @sin(theta)));
    }

    pub fn rotateVec(v: Vec3, q: Quat) Vec3 {
        const ux = q.x;
        const uy = q.y;
        const uz = q.z;
        const s = q.w;
        const vx = v.x;
        const vy = v.y;
        const vz = v.z;
        const cx = uy * vz - uz * vy;
        const cy = uz * vx - ux * vz;
        const cz = ux * vy - uy * vx;
        const c2x = uy * cz - uz * cy + s * vx;
        const c2y = uz * cx - ux * cz + s * vy;
        const c2z = ux * cy - uy * cx + s * vz;
        return Vec3{
            .x = vx + 2.0 * (uy * c2z - uz * c2y),
            .y = vy + 2.0 * (uz * c2x - ux * c2z),
            .z = vz + 2.0 * (ux * c2y - uy * c2x),
        };
    }
};

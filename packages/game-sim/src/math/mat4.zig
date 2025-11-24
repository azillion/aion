const std = @import("std");
const vec3 = @import("vec3.zig");
const quat = @import("quat.zig");

pub const Mat4 = struct {
    m: [16]f64,

    const Self = @This();

    pub fn identity() Self {
        return Self{
            .m = [_]f64{
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1,
            },
        };
    }

    pub fn fromPerspective(fovy: f64, aspect: f64, near: f64, far: f64) Self {
        std.debug.assert(aspect != 0);
        std.debug.assert(near != far);

        const t = std.math.tan(fovy * std.math.pi / 180.0 / 2.0);
        var out = Self{ .m = std.mem.zeroes([16]f64) };

        out.m[0] = 1.0 / (t * aspect);
        out.m[5] = 1.0 / t;
        out.m[10] = -(far + near) / (far - near);
        out.m[11] = -1.0;
        out.m[14] = -(2.0 * far * near) / (far - near);
        out.m[15] = 0.0;

        return out;
    }

    pub fn fromQuaternion(q: quat.Quat) Self {
        const rx = q.x;
        const ry = q.y;
        const rz = q.z;
        const rw = q.w;

        return Self{
            .m = [_]f64{ 1 - 2 * (ry * ry + rz * rz), 2 * (rx * ry - rz * rw), 2 * (rx * rz + ry * rw), 0, 2 * (rx * ry + rz * rw), 1 - 2 * (rx * rx + rz * rz), 2 * (ry * rz - rx * rw), 0, 2 * (rx * rz - ry * rw), 2 * (ry * rz + rx * rw), 1 - 2 * (rx * rx + ry * ry), 0, 0, 0, 0, 1 },
        };
    }

    pub fn lookAt(eye: vec3.Vec3, look_at: vec3.Vec3, up: vec3.Vec3) Self {
        const z_axis = look_at.sub(eye).normalize();
        const x_axis = up.cross(z_axis).normalize();
        const y_axis = z_axis.cross(x_axis);

        var out = Self{ .m = std.mem.zeroes([16]f64) };

        out.m[0] = x_axis.x;
        out.m[1] = y_axis.x;
        out.m[2] = -z_axis.x;
        out.m[3] = 0;
        out.m[4] = x_axis.y;
        out.m[5] = y_axis.y;
        out.m[6] = -z_axis.y;
        out.m[7] = 0;
        out.m[8] = x_axis.z;
        out.m[9] = y_axis.z;
        out.m[10] = -z_axis.z;
        out.m[11] = 0;
        out.m[12] = -vec3.Vec3.dot(x_axis, eye);
        out.m[13] = -vec3.Vec3.dot(y_axis, eye);
        out.m[14] = vec3.Vec3.dot(z_axis, eye);
        out.m[15] = 1;

        return out;
    }

    pub fn multiply(a: Self, b: Self) Self {
        var out = Self{ .m = std.mem.zeroes([16]f64) };

        out.m[0] = a.m[0] * b.m[0] + a.m[4] * b.m[1] + a.m[8] * b.m[2] + a.m[12] * b.m[3];
        out.m[1] = a.m[1] * b.m[0] + a.m[5] * b.m[1] + a.m[9] * b.m[2] + a.m[13] * b.m[3];
        out.m[2] = a.m[2] * b.m[0] + a.m[6] * b.m[1] + a.m[10] * b.m[2] + a.m[14] * b.m[3];
        out.m[3] = a.m[3] * b.m[0] + a.m[7] * b.m[1] + a.m[11] * b.m[2] + a.m[15] * b.m[3];

        out.m[4] = a.m[0] * b.m[4] + a.m[4] * b.m[5] + a.m[8] * b.m[6] + a.m[12] * b.m[7];
        out.m[5] = a.m[1] * b.m[4] + a.m[5] * b.m[5] + a.m[9] * b.m[6] + a.m[13] * b.m[7];
        out.m[6] = a.m[2] * b.m[4] + a.m[6] * b.m[5] + a.m[10] * b.m[6] + a.m[14] * b.m[7];
        out.m[7] = a.m[3] * b.m[4] + a.m[7] * b.m[5] + a.m[11] * b.m[6] + a.m[15] * b.m[7];

        out.m[8] = a.m[0] * b.m[8] + a.m[4] * b.m[9] + a.m[8] * b.m[10] + a.m[12] * b.m[11];
        out.m[9] = a.m[1] * b.m[8] + a.m[5] * b.m[9] + a.m[9] * b.m[10] + a.m[13] * b.m[11];
        out.m[10] = a.m[2] * b.m[8] + a.m[6] * b.m[9] + a.m[10] * b.m[10] + a.m[14] * b.m[11];
        out.m[11] = a.m[3] * b.m[8] + a.m[7] * b.m[9] + a.m[11] * b.m[10] + a.m[15] * b.m[11];

        out.m[12] = a.m[0] * b.m[12] + a.m[4] * b.m[13] + a.m[8] * b.m[14] + a.m[12] * b.m[15];
        out.m[13] = a.m[1] * b.m[12] + a.m[5] * b.m[13] + a.m[9] * b.m[14] + a.m[13] * b.m[15];
        out.m[14] = a.m[2] * b.m[12] + a.m[6] * b.m[13] + a.m[10] * b.m[14] + a.m[14] * b.m[15];
        out.m[15] = a.m[3] * b.m[12] + a.m[7] * b.m[13] + a.m[11] * b.m[14] + a.m[15] * b.m[15];

        return out;
    }

    pub fn multiplyVec3(self: Self, v: vec3.Vec3) vec3.Vec3 {
        const v4x = v.x * self.m[0] + v.y * self.m[4] + v.z * self.m[8] + self.m[12];
        const v4y = v.x * self.m[1] + v.y * self.m[5] + v.z * self.m[9] + self.m[13];
        const v4z = v.x * self.m[2] + v.y * self.m[6] + v.z * self.m[10] + self.m[14];
        const v4w = v.x * self.m[3] + v.y * self.m[7] + v.z * self.m[11] + self.m[15];
        const eps: f64 = 1e-9;
        if (@abs(v4w) < eps) {
            return vec3.Vec3{
                .x = v4x,
                .y = v4y,
                .z = v4z,
            };
        }
        const inv_w = 1.0 / v4w;

        return vec3.Vec3{
            .x = v4x * inv_w,
            .y = v4y * inv_w,
            .z = v4z * inv_w,
        };
    }

    pub fn translate(self: Self, t: vec3.Vec3) Self {
        const translation_matrix = Self{
            .m = [_]f64{
                1,   0,   0,   0,
                0,   1,   0,   0,
                0,   0,   1,   0,
                t.x, t.y, t.z, 1,
            },
        };
        return translation_matrix.multiply(self);
    }

    pub fn scale(self: Self, s: vec3.Vec3) Self {
        const scale_matrix = Self{
            .m = [_]f64{
                s.x, 0,   0,   0,
                0,   s.y, 0,   0,
                0,   0,   s.z, 0,
                0,   0,   0,   1,
            },
        };
        return scale_matrix.multiply(self);
    }

    pub fn invert(self: Self) Self {
        const a = self.m;
        var out = Self{ .m = std.mem.zeroes([16]f64) };

        out.m[0] = a[6] * a[11] * a[15] - a[6] * a[12] * a[14] - a[10] * a[7] * a[15] + a[10] * a[8] * a[14] + a[14] * a[7] * a[12] - a[14] * a[8] * a[11];
        out.m[1] = -a[2] * a[11] * a[15] + a[2] * a[12] * a[14] + a[10] * a[3] * a[15] - a[10] * a[4] * a[14] - a[14] * a[3] * a[12] + a[14] * a[4] * a[11];
        out.m[2] = a[2] * a[7] * a[15] - a[2] * a[8] * a[14] - a[6] * a[3] * a[15] + a[6] * a[4] * a[14] + a[14] * a[3] * a[8] - a[14] * a[4] * a[7];
        out.m[3] = -a[2] * a[7] * a[12] + a[2] * a[8] * a[11] + a[6] * a[3] * a[12] - a[6] * a[4] * a[11] - a[10] * a[3] * a[8] + a[10] * a[4] * a[7];

        out.m[4] = -a[5] * a[11] * a[15] + a[5] * a[12] * a[14] + a[9] * a[7] * a[15] - a[9] * a[8] * a[14] - a[13] * a[7] * a[12] + a[13] * a[8] * a[11];
        out.m[5] = a[1] * a[11] * a[15] - a[1] * a[12] * a[14] - a[9] * a[3] * a[15] + a[9] * a[4] * a[14] + a[13] * a[3] * a[12] - a[13] * a[4] * a[11];
        out.m[6] = -a[1] * a[7] * a[15] + a[1] * a[8] * a[14] + a[5] * a[3] * a[15] - a[5] * a[4] * a[14] - a[13] * a[3] * a[8] + a[13] * a[4] * a[7];
        out.m[7] = a[1] * a[7] * a[12] - a[1] * a[8] * a[11] - a[5] * a[3] * a[12] + a[5] * a[4] * a[11] + a[9] * a[3] * a[8] - a[9] * a[4] * a[7];

        out.m[8] = a[5] * a[10] * a[15] - a[5] * a[12] * a[13] - a[9] * a[6] * a[15] + a[9] * a[8] * a[13] + a[13] * a[6] * a[12] - a[13] * a[8] * a[10];
        out.m[9] = -a[1] * a[10] * a[15] + a[1] * a[12] * a[13] + a[9] * a[2] * a[15] - a[9] * a[4] * a[13] - a[13] * a[2] * a[12] + a[13] * a[4] * a[10];
        out.m[10] = a[1] * a[6] * a[15] - a[1] * a[8] * a[13] - a[5] * a[2] * a[15] + a[5] * a[4] * a[13] + a[13] * a[2] * a[8] - a[13] * a[4] * a[6];
        out.m[11] = -a[1] * a[6] * a[12] + a[1] * a[8] * a[10] + a[5] * a[2] * a[12] - a[5] * a[4] * a[10] - a[9] * a[2] * a[8] + a[9] * a[4] * a[6];

        out.m[12] = -a[5] * a[10] * a[14] + a[5] * a[11] * a[13] + a[9] * a[6] * a[14] - a[9] * a[7] * a[13] - a[12] * a[6] * a[11] + a[12] * a[7] * a[10];
        out.m[13] = a[1] * a[10] * a[14] - a[1] * a[11] * a[13] - a[9] * a[2] * a[14] + a[9] * a[3] * a[13] + a[12] * a[2] * a[11] - a[12] * a[3] * a[10];
        out.m[14] = -a[1] * a[6] * a[14] + a[1] * a[7] * a[13] + a[5] * a[2] * a[14] - a[5] * a[3] * a[13] - a[12] * a[2] * a[7] + a[12] * a[3] * a[6];
        out.m[15] = a[1] * a[6] * a[11] - a[1] * a[7] * a[10] - a[5] * a[2] * a[11] + a[5] * a[3] * a[10] + a[9] * a[2] * a[7] - a[9] * a[3] * a[6];

        const det = a[0] * out.m[0] + a[1] * out.m[4] + a[2] * out.m[8] + a[3] * out.m[12];
        const inv_det = 1.0 / det;

        for (0..16) |i| {
            out.m[i] *= inv_det;
        }

        return out;
    }

    pub fn format(
        self: Mat4,
        comptime fmt: []const u8,
        options: std.fmt.FormatOptions,
        writer: anytype,
    ) !void {
        _ = fmt;
        _ = options;
        try writer.print("[\n", .{});
        try writer.print("  {d:+.3}, {d:+.3}, {d:+.3}, {d:+.3},\n", .{ self.m[0], self.m[1], self.m[2], self.m[3] });
        try writer.print("  {d:+.3}, {d:+.3}, {d:+.3}, {d:+.3},\n", .{ self.m[4], self.m[5], self.m[6], self.m[7] });
        try writer.print("  {d:+.3}, {d:+.3}, {d:+.3}, {d:+.3},\n", .{ self.m[8], self.m[9], self.m[10], self.m[11] });
        try writer.print("  {d:+.3}, {d:+.3}, {d:+.3}, {d:+.3}\n", .{ self.m[12], self.m[13], self.m[14], self.m[15] });
        try writer.print("]", .{});
    }
};

test "Mat4 invert random affine transforms" {
    var prng = std.Random.DefaultPrng.init(0);
    var random = prng.random();

    var iteration: usize = 0;
    while (iteration < 64) : (iteration += 1) {
        const axis_raw = vec3.Vec3{
            .x = random.float(f64) * 2.0 - 1.0,
            .y = random.float(f64) * 2.0 - 1.0,
            .z = random.float(f64) * 2.0 - 1.0,
        };
        const axis = if (axis_raw.isZero()) vec3.Vec3{ .x = 0.0, .y = 0.0, .z = 1.0 } else axis_raw.normalize();
        const angle = random.float(f64) * 2.0 * std.math.pi;
        const rotation = Mat4.fromQuaternion(quat.Quat.fromAxisAngle(axis, angle));

        const scale_vec = vec3.Vec3{
            .x = random.float(f64) * 2.0 + 0.25,
            .y = random.float(f64) * 2.0 + 0.25,
            .z = random.float(f64) * 2.0 + 0.25,
        };
        const scale_mat = Mat4.identity().scale(scale_vec);

        const translation = vec3.Vec3{
            .x = random.float(f64) * 200.0 - 100.0,
            .y = random.float(f64) * 200.0 - 100.0,
            .z = random.float(f64) * 200.0 - 100.0,
        };
        const translation_mat = Mat4.identity().translate(translation);

        const transform = Mat4.multiply(translation_mat, Mat4.multiply(rotation, scale_mat));
        const inv = transform.invert();

        const ident = Mat4.multiply(transform, inv);
        const other = Mat4.multiply(inv, transform);
        inline for (0..16) |idx| {
            const expected = if (idx % 5 == 0) 1.0 else 0.0;
            try std.testing.expectApproxEqAbs(expected, ident.m[idx], 1e-9);
            try std.testing.expectApproxEqAbs(expected, other.m[idx], 1e-9);
        }
    }
}

test "Mat4 multiplyVec3 tolerates zero w" {
    var mat = Mat4.identity();
    mat.m[3] = 0.0;
    mat.m[7] = 0.0;
    mat.m[11] = 0.0;
    mat.m[15] = 0.0;

    const dir = vec3.Vec3{ .x = 2.0, .y = -3.0, .z = 5.0 };
    const result = mat.multiplyVec3(dir);

    try std.testing.expectApproxEqAbs(2.0, result.x, 1e-12);
    try std.testing.expectApproxEqAbs(-3.0, result.y, 1e-12);
    try std.testing.expectApproxEqAbs(5.0, result.z, 1e-12);
}

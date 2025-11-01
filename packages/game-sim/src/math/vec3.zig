const std = @import("std");

pub const Vec3 = struct {
    x: f64,
    y: f64,
    z: f64,

    const Self = @This();

    pub fn len(self: Vec3) f64 {
        return @sqrt(self.x * self.x + self.y * self.y + self.z * self.z);
    }

    pub fn scale(self: Vec3, s: f64) Vec3 {
        return Vec3{ .x = self.x * s, .y = self.y * s, .z = self.z * s };
    }

    pub fn add(a: Vec3, b: Vec3) Vec3 {
        return Vec3{ .x = a.x + b.x, .y = a.y + b.y, .z = a.z + b.z };
    }

    pub fn sub(a: Vec3, b: Vec3) Vec3 {
        return Vec3{ .x = a.x - b.x, .y = a.y - b.y, .z = a.z - b.z };
    }

    pub fn cross(a: Vec3, b: Vec3) Vec3 {
        return Vec3{
            .x = a.y * b.z - a.z * b.y,
            .y = a.z * b.x - a.x * b.z,
            .z = a.x * b.y - a.y * b.x,
        };
    }

    pub fn mul(a: Vec3, b: Vec3) Vec3 {
        return Vec3{
            .x = a.x * b.x,
            .y = a.y * b.y,
            .z = a.z * b.z,
        };
    }

    pub fn div(a: Vec3, b: Vec3) Vec3 {
        return Vec3{
            .x = a.x / b.x,
            .y = a.y / b.y,
            .z = a.z / b.z,
        };
    }

    pub fn normalize(self: Vec3) Vec3 {
        if (self.isZero()) {
            return ZERO;
        }
        return self.scale(1.0 / self.len());
    }

    pub fn len2(self: Vec3) f64 {
        return self.x * self.x + self.y * self.y + self.z * self.z;
    }

    pub fn negate(self: Vec3) Vec3 {
        return Vec3{
            .x = -self.x,
            .y = -self.y,
            .z = -self.z,
        };
    }

    pub fn dot(a: Vec3, b: Vec3) f64 {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    pub fn dist(a: Vec3, b: Vec3) f64 {
        const dx = a.x - b.x;
        const dy = a.y - b.y;
        const dz = a.z - b.z;
        return @sqrt(dx * dx + dy * dy + dz * dz);
    }

    pub fn dist2(a: Vec3, b: Vec3) f64 {
        const dx = a.x - b.x;
        const dy = a.y - b.y;
        const dz = a.z - b.z;
        return dx * dx + dy * dy + dz * dz;
    }

    pub fn rotate(self: Vec3, phi: f64, axis: Vec3) Vec3 {
        const u = axis.normalize();
        const c = @cos(phi);
        const s = @sin(phi);

        // Calculate generalized rotation matrix components
        const m1_x = c + u.x * u.x * (1 - c);
        const m1_y = u.x * u.y * (1 - c) - u.z * s;
        const m1_z = u.x * u.z * (1 - c) + u.y * s;

        const m2_x = u.y * u.x * (1 - c) + u.z * s;
        const m2_y = c + u.y * u.y * (1 - c);
        const m2_z = u.y * u.z * (1 - c) - u.x * s;

        const m3_x = u.z * u.x * (1 - c) - u.y * s;
        const m3_y = u.z * u.y * (1 - c) + u.x * s;
        const m3_z = c + u.z * u.z * (1 - c);

        return Vec3{
            .x = self.x * m1_x + self.y * m1_y + self.z * m1_z,
            .y = self.x * m2_x + self.y * m2_y + self.z * m2_z,
            .z = self.x * m3_x + self.y * m3_y + self.z * m3_z,
        };
    }

    pub fn lerp(a: Vec3, b: Vec3, s: f64) Vec3 {
        return a.add(b.sub(a).scale(s));
    }

    pub fn componentMin(a: Vec3, b: Vec3) Vec3 {
        return Vec3{
            .x = @min(a.x, b.x),
            .y = @min(a.y, b.y),
            .z = @min(a.z, b.z),
        };
    }

    pub fn componentMax(a: Vec3, b: Vec3) Vec3 {
        return Vec3{
            .x = @max(a.x, b.x),
            .y = @max(a.y, b.y),
            .z = @max(a.z, b.z),
        };
    }

    pub fn isZero(self: Vec3) bool {
        return self.x == 0 and self.y == 0 and self.z == 0;
    }

    pub fn format(
        self: Vec3,
        comptime fmt: []const u8,
        options: std.fmt.FormatOptions,
        writer: anytype,
    ) !void {
        _ = fmt;
        _ = options;
        try writer.print("({d:+.3},{d:+.3},{d:+.3})", .{ self.x, self.y, self.z });
    }

    pub const UNIT_X = Vec3{ .x = 1, .y = 0, .z = 0 };
    pub const UNIT_Y = Vec3{ .x = 0, .y = 1, .z = 0 };
    pub const UNIT_Z = Vec3{ .x = 0, .y = 0, .z = 1 };
    pub const ZERO = Vec3{ .x = 0, .y = 0, .z = 0 };
};

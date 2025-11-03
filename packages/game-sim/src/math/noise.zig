const std = @import("std");
const Vec3 = @import("vec3.zig").Vec3;

// 3D Simplex Noise implementation (Stefan Gustavson style)
// Returns values roughly in [-1, 1]. Deterministic for a given seed.
pub const SimplexNoise = struct {
    perm: [512]u8,

    pub fn init(seed: u64) SimplexNoise {
        var self: SimplexNoise = undefined;
        var p: [256]u8 = undefined;
        var i: usize = 0;
        while (i < 256) : (i += 1) p[i] = @intCast(i);

        var prng = std.Random.DefaultPrng.init(seed);
        const random = prng.random();
        random.shuffle(u8, &p);

        i = 0;
        while (i < 256) : (i += 1) {
            self.perm[i] = p[i];
            self.perm[i + 256] = p[i];
        }
        return self;
    }

    inline fn dot(g: Vec3, x: f64, y: f64, z: f64) f64 {
        return g.x * x + g.y * y + g.z * z;
    }

    pub fn noise(self: *const SimplexNoise, xin: f64, yin: f64, zin: f64) f64 {
        // Skewing and unskewing factors for 3D
        const F3 = 1.0 / 3.0;
        const G3 = 1.0 / 6.0;

        const s = (xin + yin + zin) * F3;
        const i_f = @floor(xin + s);
        const j_f = @floor(yin + s);
        const k_f = @floor(zin + s);

        const i: i32 = @intFromFloat(i_f);
        const j: i32 = @intFromFloat(j_f);
        const k: i32 = @intFromFloat(k_f);

        const t = (@as(f64, @floatFromInt(i + j + k))) * G3;
        const X0 = @as(f64, @floatFromInt(i)) - t;
        const Y0 = @as(f64, @floatFromInt(j)) - t;
        const Z0 = @as(f64, @floatFromInt(k)) - t;
        const x0 = xin - X0;
        const y0 = yin - Y0;
        const z0 = zin - Z0;

        var i1s: i32 = 0;
        var j1s: i32 = 0;
        var k1s: i32 = 0;
        var i2s: i32 = 0;
        var j2s: i32 = 0;
        var k2s: i32 = 0;

        if (x0 >= y0) {
            if (y0 >= z0) {
                i1s = 1;
                j1s = 0;
                k1s = 0;
                i2s = 1;
                j2s = 1;
                k2s = 0;
            } else if (x0 >= z0) {
                i1s = 1;
                j1s = 0;
                k1s = 0;
                i2s = 1;
                j2s = 0;
                k2s = 1;
            } else {
                i1s = 0;
                j1s = 0;
                k1s = 1;
                i2s = 1;
                j2s = 0;
                k2s = 1;
            }
        } else {
            if (y0 < z0) {
                i1s = 0;
                j1s = 0;
                k1s = 1;
                i2s = 0;
                j2s = 1;
                k2s = 1;
            } else if (x0 < z0) {
                i1s = 0;
                j1s = 1;
                k1s = 0;
                i2s = 0;
                j2s = 1;
                k2s = 1;
            } else {
                i1s = 0;
                j1s = 1;
                k1s = 0;
                i2s = 1;
                j2s = 1;
                k2s = 0;
            }
        }

        const x1 = x0 - @as(f64, @floatFromInt(i1s)) + G3;
        const y1 = y0 - @as(f64, @floatFromInt(j1s)) + G3;
        const z1 = z0 - @as(f64, @floatFromInt(k1s)) + G3;
        const x2 = x0 - @as(f64, @floatFromInt(i2s)) + 2.0 * G3;
        const y2 = y0 - @as(f64, @floatFromInt(j2s)) + 2.0 * G3;
        const z2 = z0 - @as(f64, @floatFromInt(k2s)) + 2.0 * G3;
        const x3 = x0 - 1.0 + 3.0 * G3;
        const y3 = y0 - 1.0 + 3.0 * G3;
        const z3 = z0 - 1.0 + 3.0 * G3;

        const grad3 = [_]Vec3{
            .{ .x = 1, .y = 1, .z = 0 }, .{ .x = -1, .y = 1, .z = 0 }, .{ .x = 1, .y = -1, .z = 0 }, .{ .x = -1, .y = -1, .z = 0 },
            .{ .x = 1, .y = 0, .z = 1 }, .{ .x = -1, .y = 0, .z = 1 }, .{ .x = 1, .y = 0, .z = -1 }, .{ .x = -1, .y = 0, .z = -1 },
            .{ .x = 0, .y = 1, .z = 1 }, .{ .x = 0, .y = -1, .z = 1 }, .{ .x = 0, .y = 1, .z = -1 }, .{ .x = 0, .y = -1, .z = -1 },
        };

        const ii: usize = @intCast(@as(i32, i) & 255);
        const jj: usize = @intCast(@as(i32, j) & 255);
        const kk: usize = @intCast(@as(i32, k) & 255);

        // Helper indices with wrapping
        const ii_i1 = (ii + @as(usize, @intCast(i1s))) & 255;
        const jj_j1 = (jj + @as(usize, @intCast(j1s))) & 255;
        const kk_k1 = (kk + @as(usize, @intCast(k1s))) & 255;
        const ii_i2 = (ii + @as(usize, @intCast(i2s))) & 255;
        const jj_j2 = (jj + @as(usize, @intCast(j2s))) & 255;
        const kk_k2 = (kk + @as(usize, @intCast(k2s))) & 255;
        const ii_1 = (ii + 1) & 255;
        const jj_1 = (jj + 1) & 255;
        const kk_1 = (kk + 1) & 255;

        const gi0 = self.perm[ii + @as(usize, self.perm[jj + @as(usize, self.perm[kk])])] % 12;
        const gi1 = self.perm[ii_i1 + @as(usize, self.perm[jj_j1 + @as(usize, self.perm[kk_k1])])] % 12;
        const gi2 = self.perm[ii_i2 + @as(usize, self.perm[jj_j2 + @as(usize, self.perm[kk_k2])])] % 12;
        const gi3 = self.perm[ii_1 + @as(usize, self.perm[jj_1 + @as(usize, self.perm[kk_1])])] % 12;

        var n0: f64 = 0.0;
        var n1: f64 = 0.0;
        var n2: f64 = 0.0;
        var n3: f64 = 0.0;

        var t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
        if (t0 > 0) {
            t0 *= t0;
            n0 = t0 * t0 * dot(grad3[gi0], x0, y0, z0);
        }

        var t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
        if (t1 > 0) {
            t1 *= t1;
            n1 = t1 * t1 * dot(grad3[gi1], x1, y1, z1);
        }

        var t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
        if (t2 > 0) {
            t2 *= t2;
            n2 = t2 * t2 * dot(grad3[gi2], x2, y2, z2);
        }

        var t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
        if (t3 > 0) {
            t3 *= t3;
            n3 = t3 * t3 * dot(grad3[gi3], x3, y3, z3);
        }

        // Scale to roughly [-1, 1]
        return 32.0 * (n0 + n1 + n2 + n3);
    }
};

const expect = std.testing.expect;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;

test "simplex: determinism and rough bounds" {
    var n1 = SimplexNoise.init(42);
    var n2 = SimplexNoise.init(42);
    var n3 = SimplexNoise.init(43);

    const a = n1.noise(0.1, 0.2, 0.3);
    const b = n2.noise(0.1, 0.2, 0.3);
    const c = n3.noise(0.1, 0.2, 0.3);

    try expectApproxEqAbs(a, b, 1e-12);
    // Different seeds should usually differ at the same sample point
    try expect(a != c);
    // Reasonable bounds for Gustavson 3D simplex (scaled by 32)
    try expect(a <= 1.25 and a >= -1.25);
}

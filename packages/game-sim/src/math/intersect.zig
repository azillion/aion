const std = @import("std");
const Vec3 = @import("vec3.zig").Vec3;
const Bound3 = @import("bound3.zig").Bound3;

pub const Ray = struct {
    origin: Vec3,
    direction: Vec3,
};

pub const Sphere = struct {
    center: Vec3,
    radius: f64,
};

pub fn aabbAabb(a: Bound3, b: Bound3) bool {
    return a.min.x <= b.max.x and
        a.max.x >= b.min.x and
        a.min.y <= b.max.y and
        a.max.y >= b.min.y and
        a.min.z <= b.max.z and
        a.max.z >= b.min.z;
}

pub fn aabbSphere(aabb: Bound3, sphere: Sphere) bool {
    var dist_sq = sphere.radius * sphere.radius;

    // Check each axis
    if (sphere.center.x < aabb.min.x) {
        const diff = sphere.center.x - aabb.min.x;
        dist_sq -= diff * diff;
    } else if (sphere.center.x > aabb.max.x) {
        const diff = sphere.center.x - aabb.max.x;
        dist_sq -= diff * diff;
    }

    if (sphere.center.y < aabb.min.y) {
        const diff = sphere.center.y - aabb.min.y;
        dist_sq -= diff * diff;
    } else if (sphere.center.y > aabb.max.y) {
        const diff = sphere.center.y - aabb.max.y;
        dist_sq -= diff * diff;
    }

    if (sphere.center.z < aabb.min.z) {
        const diff = sphere.center.z - aabb.min.z;
        dist_sq -= diff * diff;
    } else if (sphere.center.z > aabb.max.z) {
        const diff = sphere.center.z - aabb.max.z;
        dist_sq -= diff * diff;
    }

    return dist_sq >= 0;
}

inline fn axisIntersect(
    origin: f64,
    dir: f64,
    min_val: f64,
    max_val: f64,
    tmin_ref: *f64,
    tmax_ref: *f64,
    eps: f64,
) bool {
    if (@abs(dir) < eps) {
        return !(origin < min_val or origin > max_val);
    }
    const inv = 1.0 / dir;
    const t1 = (min_val - origin) * inv;
    const t2 = (max_val - origin) * inv;
    tmin_ref.* = @max(tmin_ref.*, @min(t1, t2));
    tmax_ref.* = @min(tmax_ref.*, @max(t1, t2));
    return tmin_ref.* <= tmax_ref.*;
}

pub fn rayAabb(ray: Ray, aabb: Bound3) ?f64 {
    const eps: f64 = 1e-12;
    const d = ray.direction;
    var tmin = -std.math.inf(f64);
    var tmax = std.math.inf(f64);

    if (!axisIntersect(ray.origin.x, d.x, aabb.min.x, aabb.max.x, &tmin, &tmax, eps)) return null;
    if (!axisIntersect(ray.origin.y, d.y, aabb.min.y, aabb.max.y, &tmin, &tmax, eps)) return null;
    if (!axisIntersect(ray.origin.z, d.z, aabb.min.z, aabb.max.z, &tmin, &tmax, eps)) return null;

    if (tmax < 0) return null;
    if (tmin > tmax) return null;
    return tmin;
}

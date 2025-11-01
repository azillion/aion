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

    return dist_sq > 0;
}

pub fn rayAabb(ray: Ray, aabb: Bound3) ?f64 {
    const dir = ray.direction.normalize();
    const dirfrac = Vec3.init(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);

    const t1 = (aabb.min.x - ray.origin.x) * dirfrac.x;
    const t2 = (aabb.max.x - ray.origin.x) * dirfrac.x;
    const t3 = (aabb.min.y - ray.origin.y) * dirfrac.y;
    const t4 = (aabb.max.y - ray.origin.y) * dirfrac.y;
    const t5 = (aabb.min.z - ray.origin.z) * dirfrac.z;
    const t6 = (aabb.max.z - ray.origin.z) * dirfrac.z;

    const tmin = @max(@max(@min(t1, t2), @min(t3, t4)), @min(t5, t6));
    const tmax = @min(@min(@max(t1, t2), @max(t3, t4)), @max(t5, t6));

    // ray is intersecting AABB, but whole AABB is behind us
    if (tmax < 0) {
        return null;
    }

    // ray does not intersect AABB
    if (tmin > tmax) {
        return null;
    }

    return tmin;
}

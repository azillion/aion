const std = @import("std");
const Vec3 = @import("vec3.zig").Vec3;

pub const Bound3 = struct {
    min: Vec3,
    max: Vec3,

    pub fn init(min: Vec3, max: Vec3) Bound3 {
        return Bound3{
            .min = min,
            .max = max,
        };
    }

    pub fn center(self: Bound3) Vec3 {
        return self.min.add(self.max).scale(0.5);
    }

    pub fn size(self: Bound3) Vec3 {
        return self.max.sub(self.min);
    }

    pub fn extend(self: Bound3, point: Vec3) Bound3 {
        return Bound3.init(
            self.min.componentMin(point),
            self.max.componentMax(point),
        );
    }

    pub fn offset(self: Bound3, offset_vec: Vec3) Bound3 {
        return Bound3.init(
            self.min.add(offset_vec),
            self.max.add(offset_vec),
        );
    }

    pub fn contains(self: Bound3, point: Vec3) bool {
        return self.min.x <= point.x and self.min.y <= point.y and self.min.z <= point.z and
            self.max.x >= point.x and self.max.y >= point.y and self.max.z >= point.z;
    }

    pub const zero = Bound3.init(Vec3.zero, Vec3.zero);
};

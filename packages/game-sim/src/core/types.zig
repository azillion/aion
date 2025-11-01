const std = @import("std");
const math = @import("math");

pub const Vec3 = math.Vec3;
pub const Quat = math.Quat;
pub const Atmosphere = struct {
    N2: f64,
    O2: f64,
};
pub const TerrainParams = struct {
    radius: f64,
    seaLevel: f64,
    maxHeight: f64,
    noiseSeed: f64,
    atmosphere: ?Atmosphere = null,
};

pub const Body = struct {
    const Self = @This();
    id: []const u8,
    name: []const u8,
    position: Vec3,
    velocity: Vec3,
    radius: f64,
    mass: f64,
    albedo: Vec3,
    emissive: ?Vec3 = null,
    terrain: ?TerrainParams = null,
};

pub const Ship = struct {
    const Self = @This();
    body: Body,

    // Ship-specific fields
    orientation: Quat,
    angularVelocity: Vec3,
    thrust: Vec3,
};

pub const SystemState = struct {
    const Self = @This();
    timestamp: f64,
    bodies: []Body,
    ships: []Ship,
};

pub const InputState = struct {
    const Self = @This();
    keys_mask: u32,
};

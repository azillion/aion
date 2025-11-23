const prebaked = @import("prebaked.zig");

pub const TileId = prebaked.TileId;
pub const UNSET_NEIGHBOR = prebaked.UNSET_NEIGHBOR;

pub const magic_u32: u32 =
    (@as(u32, 'P')) |
    (@as(u32, 'L') << 8) |
    (@as(u32, 'N') << 16) |
    (@as(u32, 'T') << 24);

pub const version: u16 = 2;
pub const flag_is_pentagon: u8 = 1;

/// Header describing a prebaked planet blob on disk (little-endian).
pub const PlanetHeader = extern struct {
    magic: u32 = magic_u32,
    version: u16 = version,
    _pad: u16 = 0,
    R: u32,
    tile_count: u32,
    vertex_count: u32,
    index_count: u32,
    edge_index_count: u32,
};

/// Packed representation of a `Tile` on disk.
pub const DiskTile = extern struct {
    pos_sphere: [3]f32,
    face: u8,
    _pad0: [3]u8 = .{0} ** 3,
    axial_q: i32,
    axial_r: i32,
    neighbors: [6]TileId,
    flags: u8,
    _pad1: [3]u8 = .{0} ** 3,
};

pub fn tileToDisk(t: prebaked.Tile) DiskTile {
    return .{
        .pos_sphere = t.pos_sphere,
        .face = t.face,
        .axial_q = t.axial.q,
        .axial_r = t.axial.r,
        .neighbors = t.neighbors,
        .flags = if (t.is_pentagon) flag_is_pentagon else 0,
    };
}

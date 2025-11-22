pub const TileId = u32;

pub const UNSET_NEIGHBOR: TileId = 0xFFFFFFFF;

pub const Tile = struct {
    pos_sphere: [3]f32, // xyz on unit sphere
    face: u8, // original face (debug/viz)
    axial: struct { q: i32, r: i32 }, // original coords (debug/viz)
    neighbors: [6]TileId, // Neighbor IDs. 0xFFFFFFFF if missing (pentagon)
    is_pentagon: bool,
};

pub const PlanetGrid = struct {
    tiles: []const Tile,
};

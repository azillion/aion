const std = @import("std");
const prebaked = @import("prebaked.zig");
const format = @import("prebaked_format.zig");
const TileId = prebaked.TileId;

pub const LoadError = error{
    BadPlanetFile,
    BadMagic,
    BadVersion,
    Truncated,
};

pub const LoadedPlanet = struct {
    header: format.PlanetHeader,
    tiles: []prebaked.Tile,

    pub fn asGrid(self: *const LoadedPlanet) prebaked.PlanetGrid {
        return .{ .tiles = self.tiles };
    }
};

pub fn loadPlanetGrid(alloc: std.mem.Allocator, bytes: []const u8) !LoadedPlanet {
    if (bytes.len < @sizeOf(format.PlanetHeader)) return LoadError.BadPlanetFile;

    const header = std.mem.bytesAsValue(format.PlanetHeader, bytes[0..@sizeOf(format.PlanetHeader)]).*;
    if (header.magic != format.magic_u32) return LoadError.BadMagic;
    if (header.version != format.version) return LoadError.BadVersion;

    const tile_count: usize = @intCast(header.tile_count);
    const expected_len = @sizeOf(format.PlanetHeader) + tile_count * @sizeOf(format.DiskTile);
    if (bytes.len < expected_len) return LoadError.Truncated;

    var tiles = try alloc.alloc(prebaked.Tile, tile_count);
    errdefer alloc.free(tiles);

    const tile_bytes = bytes[@sizeOf(format.PlanetHeader)..expected_len];
    var i: usize = 0;
    while (i < tile_count) : (i += 1) {
        const offset = i * @sizeOf(format.DiskTile);
        const disk_slice = tile_bytes[offset .. offset + @sizeOf(format.DiskTile)];
        const disk = std.mem.bytesAsValue(format.DiskTile, disk_slice).*;
        tiles[i] = .{
            .pos_sphere = disk.pos_sphere,
            .face = disk.face,
            .axial = .{ .q = disk.axial_q, .r = disk.axial_r },
            .neighbors = disk.neighbors,
            .is_pentagon = (disk.flags & format.flag_is_pentagon) != 0,
        };
    }

    return .{
        .header = header,
        .tiles = tiles,
    };
}

test "loadPlanetGrid roundtrip" {
    const alloc = std.testing.allocator;

    const header = format.PlanetHeader{
        .R = 4,
        .tile_count = 1,
    };

    const disk_tile = format.DiskTile{
        .pos_sphere = .{ 0.0, 1.0, 0.0 },
        .face = 2,
        .axial_q = -1,
        .axial_r = 3,
        .neighbors = .{ 1, 2, 3, 4, 5, 6 },
        .flags = format.flag_is_pentagon,
    };

    var blob: [@sizeOf(format.PlanetHeader) + @sizeOf(format.DiskTile)]u8 = undefined;
    const header_bytes = std.mem.asBytes(&header);
    const tile_bytes = std.mem.asBytes(&disk_tile);
    std.mem.copyForwards(u8, blob[0..header_bytes.len], header_bytes);
    std.mem.copyForwards(u8, blob[header_bytes.len..], tile_bytes);

    const asset = try loadPlanetGrid(alloc, &blob);
    defer alloc.free(asset.tiles);

    try std.testing.expectEqual(@as(usize, 1), asset.tiles.len);
    const tile = asset.tiles[0];
    try std.testing.expect(tile.is_pentagon);
    try std.testing.expectEqual(@as(u8, 2), tile.face);
    try std.testing.expectEqual(@as(i32, -1), tile.axial.q);
    try std.testing.expectEqual(@as(TileId, 4), tile.neighbors[3]);
}

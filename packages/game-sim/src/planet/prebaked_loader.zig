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
    vertices: []f32,
    elevations: []f32,
    indices: []u32,
    edges: []u32,

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
    const vertex_count: usize = @intCast(header.vertex_count);
    const index_count: usize = @intCast(header.index_count);
    const edge_count: usize = @intCast(header.edge_index_count);
    const vertex_f32_count = vertex_count * 3;
    const elevation_f32_count = vertex_count;

    const tile_bytes_len = tile_count * @sizeOf(format.DiskTile);
    const vertex_bytes_len = vertex_f32_count * @sizeOf(f32);
    const elevation_bytes_len = elevation_f32_count * @sizeOf(f32);
    const index_bytes_len = index_count * @sizeOf(u32);
    const edge_bytes_len = edge_count * @sizeOf(u32);

    var offset: usize = @sizeOf(format.PlanetHeader);
    const total_len = offset + tile_bytes_len + vertex_bytes_len + elevation_bytes_len + index_bytes_len + edge_bytes_len;
    if (bytes.len < total_len) return LoadError.Truncated;

    var tiles = try alloc.alloc(prebaked.Tile, tile_count);
    errdefer alloc.free(tiles);

    const vertices = try alloc.alloc(f32, vertex_f32_count);
    errdefer alloc.free(vertices);
    const elevations = try alloc.alloc(f32, elevation_f32_count);
    errdefer alloc.free(elevations);
    const indices = try alloc.alloc(u32, index_count);
    errdefer alloc.free(indices);
    const edges = try alloc.alloc(u32, edge_count);
    errdefer alloc.free(edges);

    const tile_region = bytes[offset .. offset + tile_bytes_len];
    offset += tile_bytes_len;
    var i: usize = 0;
    while (i < tile_count) : (i += 1) {
        const base = i * @sizeOf(format.DiskTile);
        const disk_slice = tile_region[base .. base + @sizeOf(format.DiskTile)];
        const disk = std.mem.bytesAsValue(format.DiskTile, disk_slice).*;
        tiles[i] = .{
            .pos_sphere = disk.pos_sphere,
            .face = disk.face,
            .axial = .{ .q = disk.axial_q, .r = disk.axial_r },
            .neighbors = disk.neighbors,
            .is_pentagon = (disk.flags & format.flag_is_pentagon) != 0,
        };
    }

    const vertex_region = bytes[offset .. offset + vertex_bytes_len];
    std.mem.copyForwards(u8, std.mem.sliceAsBytes(vertices), vertex_region);
    offset += vertex_bytes_len;

    const elevation_region = bytes[offset .. offset + elevation_bytes_len];
    std.mem.copyForwards(u8, std.mem.sliceAsBytes(elevations), elevation_region);
    offset += elevation_bytes_len;

    const index_region = bytes[offset .. offset + index_bytes_len];
    std.mem.copyForwards(u8, std.mem.sliceAsBytes(indices), index_region);
    offset += index_bytes_len;

    const edge_region = bytes[offset .. offset + edge_bytes_len];
    std.mem.copyForwards(u8, std.mem.sliceAsBytes(edges), edge_region);

    return .{
        .header = header,
        .tiles = tiles,
        .vertices = vertices,
        .elevations = elevations,
        .indices = indices,
        .edges = edges,
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
    defer {
        alloc.free(asset.tiles);
        alloc.free(asset.vertices);
        alloc.free(asset.elevations);
        alloc.free(asset.indices);
        alloc.free(asset.edges);
    }

    try std.testing.expectEqual(@as(usize, 1), asset.tiles.len);
    const tile = asset.tiles[0];
    try std.testing.expect(tile.is_pentagon);
    try std.testing.expectEqual(@as(u8, 2), tile.face);
    try std.testing.expectEqual(@as(i32, -1), tile.axial.q);
    try std.testing.expectEqual(@as(TileId, 4), tile.neighbors[3]);
    try std.testing.expectEqual(@as(usize, 0), asset.vertices.len);
    try std.testing.expectEqual(@as(usize, 0), asset.indices.len);
    try std.testing.expectEqual(@as(usize, 0), asset.edges.len);
}

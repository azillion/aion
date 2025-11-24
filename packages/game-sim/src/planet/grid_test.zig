const std = @import("std");
const grid_mod = @import("grid.zig");
const builder = @import("prebaked_builder.zig");
const prebaked = @import("prebaked.zig");

const Grid = grid_mod.Grid;
const Tile = prebaked.Tile;

test "hashCoord roundtrip" {
    const rand_count: usize = 2000;
    var prng = std.Random.DefaultPrng.init(0xC0FFEE);
    var random = prng.random();

    var i: usize = 0;
    while (i < rand_count) : (i += 1) {
        const q: isize = @intCast(random.intRangeLessThan(i32, -500_000, 500_000));
        const r: isize = @intCast(random.intRangeLessThan(i32, -500_000, 500_000));
        const face: usize = @intCast(random.intRangeLessThan(u8, 0, 20));

        const hash = Grid.hashCoord(q, r, face);
        const roundtrip = Grid.unhashCoord(hash);

        try std.testing.expectEqual(q, roundtrip.q);
        try std.testing.expectEqual(r, roundtrip.r);
        try std.testing.expectEqual(face, roundtrip.face);
    }
}

test "grid neighbor rings are symmetric, unique, and ordered" {
    var grid = try Grid.init(std.testing.allocator, 3);
    defer grid.deinit();

    const baked = try builder.build_prebaked(&grid, std.testing.allocator);
    defer {
        std.testing.allocator.free(baked.tiles);
        std.testing.allocator.free(baked.vertices);
        std.testing.allocator.free(baked.elevations);
        std.testing.allocator.free(baked.indices);
        std.testing.allocator.free(baked.edges);
    }

    for (baked.tiles, 0..) |tile, idx| {
        const neighbor_count = countNeighbors(tile);
        if (tile.is_pentagon) {
            try std.testing.expectEqual(@as(usize, 5), neighbor_count);
        } else {
            try std.testing.expectEqual(@as(usize, 6), neighbor_count);
        }

        try ensureUniqueNeighbors(&tile, neighbor_count);
        try ensureSymmetricNeighbors(idx, tile, baked.tiles);
        try ensureRingOrdering(idx, tile, baked.tiles);
    }
}

fn countNeighbors(tile: Tile) usize {
    var count: usize = 0;
    for (tile.neighbors) |nb| {
        if (nb != prebaked.UNSET_NEIGHBOR) count += 1;
    }
    return count;
}

fn ensureUniqueNeighbors(tile: *const Tile, neighbor_count: usize) !void {
    var i: usize = 0;
    while (i < neighbor_count) : (i += 1) {
        const a = tile.neighbors[i];
        try std.testing.expect(a != prebaked.UNSET_NEIGHBOR);
        var j = i + 1;
        while (j < neighbor_count) : (j += 1) {
            try std.testing.expect(tile.neighbors[j] != a);
        }
    }
}

fn ensureSymmetricNeighbors(tile_idx: usize, tile: Tile, tiles: []const Tile) !void {
    for (tile.neighbors) |nb| {
        if (nb == prebaked.UNSET_NEIGHBOR) continue;
        const nb_idx: usize = @intCast(nb);
        try std.testing.expect(nb_idx < tiles.len);
        const neighbor = tiles[nb_idx];
        var seen_back = false;
        for (neighbor.neighbors) |rev| {
            if (rev == prebaked.UNSET_NEIGHBOR) continue;
            if (rev == tile_idx) {
                seen_back = true;
                break;
            }
        }
        try std.testing.expect(seen_back);
    }
}

fn ensureRingOrdering(tile_idx: usize, tile: Tile, tiles: []const Tile) !void {
    _ = tile_idx;
    var neighbor_indices: [6]usize = undefined;
    var count: usize = 0;
    for (tile.neighbors) |nb| {
        if (nb == prebaked.UNSET_NEIGHBOR) continue;
        neighbor_indices[count] = @intCast(nb);
        count += 1;
    }
    if (count <= 2) return;

    const basis = makeBasis(tile.pos_sphere);
    var angles: [6]f64 = undefined;
    var i: usize = 0;
    while (i < count) : (i += 1) {
        const nb = tiles[neighbor_indices[i]].pos_sphere;
        const rel = [3]f64{
            @floatCast(nb[0] - tile.pos_sphere[0]),
            @floatCast(nb[1] - tile.pos_sphere[1]),
            @floatCast(nb[2] - tile.pos_sphere[2]),
        };
        const x = rel[0] * basis.u[0] + rel[1] * basis.u[1] + rel[2] * basis.u[2];
        const y = rel[0] * basis.v[0] + rel[1] * basis.v[1] + rel[2] * basis.v[2];
        angles[i] = std.math.atan2(y, x);
    }

    var prev = angles[0];
    var idx: usize = 1;
    while (idx < count) : (idx += 1) {
        try std.testing.expect(angles[idx] + 1e-9 >= prev);
        if (angles[idx] > prev) {
            prev = angles[idx];
        }
    }
}

fn normalize(v: [3]f64) [3]f64 {
    const len = std.math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (len == 0) return .{ 0.0, 0.0, 0.0 };
    return .{ v[0] / len, v[1] / len, v[2] / len };
}

fn cross(a: [3]f64, b: [3]f64) [3]f64 {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

fn makeBasis(center: [3]f32) struct { u: [3]f64, v: [3]f64 } {
    const c = normalize(.{
        @as(f64, center[0]),
        @as(f64, center[1]),
        @as(f64, center[2]),
    });
    var up: [3]f64 = .{ 0.0, 0.0, 1.0 };
    var u = cross(c, up);
    if (std.math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) < 1e-6) {
        up = .{ 0.0, 1.0, 0.0 };
        u = cross(c, up);
    }
    u = normalize(u);
    const v = normalize(cross(c, u));
    return .{ .u = u, .v = v };
}

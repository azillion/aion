const std = @import("std");
const Grid = @import("grid.zig").Grid;
const prebaked = @import("prebaked.zig");
const Tile = prebaked.Tile;
const TileId = prebaked.TileId;

const DIRS = [_][2]i32{
    .{ 1, 0 },
    .{ 1, -1 },
    .{ 0, -1 },
    .{ -1, 0 },
    .{ -1, 1 },
    .{ 0, 1 },
};

const Key = struct {
    face: u8,
    q: i32,
    r: i32,
};

const CoordMap = std.AutoHashMap(Key, TileId);

fn makeKey(face: u8, q: i32, r: i32) Key {
    return .{ .face = face, .q = q, .r = r };
}

pub fn build_prebaked(grid: *Grid, alloc: std.mem.Allocator) ![]Tile {
    // Analytic per-face enumeration; no BFS or canonicalization.
    const N: isize = @intCast(grid.size);

    var tiles_list = std.ArrayListUnmanaged(Tile){};
    defer tiles_list.deinit(alloc);

    // Iterate each face independently
    for (0..20) |face| {
        var q: isize = -N;
        while (q <= N) : (q += 1) {
            var r: isize = -N;
            while (r <= N) : (r += 1) {
                // Face-local barycentric (scaled) for the axial (q,r)
                const uvw = Grid.toBary_public(N, q, r);

                // Inside the hex window of this face?
                if (uvw[0] < 0 or uvw[1] < 0 or uvw[2] < 0) continue;
                // Optional tighter guard; keep for symmetry with other helpers.
                if (uvw[0] > 2 * N or uvw[1] > 2 * N or uvw[2] > 2 * N) continue;

                const pos = grid.faceCenterPosition(q, r, face);
                const is_penta = grid.isPentaAt(face, q, r);

                try tiles_list.append(alloc, .{
                    .pos_sphere = .{ @floatCast(pos.x), @floatCast(pos.y), @floatCast(pos.z) },
                    .face = @intCast(face),
                    .axial = .{ .q = @intCast(q), .r = @intCast(r) },
                    .neighbors = .{prebaked.UNSET_NEIGHBOR} ** 6,
                    .is_pentagon = is_penta,
                });
            }
        }
    }

    const tiles = try tiles_list.toOwnedSlice(alloc);

    // Sanity: count per-face tiles and pentagons.
    var per_face = [_]usize{0} ** 20;
    var pentas: usize = 0;
    for (tiles) |t| {
        per_face[t.face] += 1;
        if (t.is_pentagon) pentas += 1;
    }
    std.debug.print("PREBAKE: per-face counts {any}\n", .{per_face});
    std.debug.print("PREBAKE: pentagons {d}\n", .{pentas});

    // Build coord -> TileId index.
    var map = CoordMap.init(alloc);
    defer map.deinit();
    for (tiles, 0..) |t, idx| {
        const key = makeKey(t.face, t.axial.q, t.axial.r);
        try map.put(key, @intCast(idx));
    }

    // Fill neighbor references.
    var ti: usize = 0;
    while (ti < tiles.len) : (ti += 1) {
        var t = &tiles[ti];
        for (DIRS, 0..) |dir, dir_idx| {
            const nq = t.axial.q + dir[0];
            const nr = t.axial.r + dir[1];
            const key = makeKey(t.face, nq, nr);
            if (map.get(key)) |nb_id| {
                t.neighbors[dir_idx] = nb_id;
                continue;
            }

            // Cross-face neighbors are deferred; leave as UNSET for now.
        }
    }

    return tiles;
}

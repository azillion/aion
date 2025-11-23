const std = @import("std");
const grid_mod = @import("grid.zig");
const Grid = grid_mod.Grid;
const icosa_face_vertices = grid_mod.icosa_face_vertices;
const prebaked = @import("prebaked.zig");
const Tile = prebaked.Tile;
const TileId = prebaked.TileId;
const UNSET_NEIGHBOR = prebaked.UNSET_NEIGHBOR;

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
const CornerMap = std.AutoHashMap([2]u8, TileId);

const BuilderTile = struct {
    tile: Tile,
    neighbor_count: u8,
};

const EdgeEntry = struct {
    tile_index: usize,
    key: isize,
};
const EdgeEntryGrid = [20][3]std.ArrayListUnmanaged(EdgeEntry);
const CornerEntry = struct {
    face: usize,
    tile_index: usize,
    vertex: u8,
};
const CornerLists = [12]std.ArrayListUnmanaged(CornerEntry);

const FaceEdgeAdj = struct {
    nf: usize,
    ne: usize,
    rev: bool,
};
const FaceEdgeAdjTable = [20][3]FaceEdgeAdj;
const face_edge_adj: FaceEdgeAdjTable = buildFaceEdgeAdjacency();

const BuiltPlanet = struct {
    tiles: []Tile,
    vertices: []f32,
    elevations: []f32,
    indices: []u32,
    edges: []u32,
};

fn makeKey(face: u8, q: i32, r: i32) Key {
    return .{ .face = face, .q = q, .r = r };
}

fn buildFaceEdgeAdjacency() FaceEdgeAdjTable {
    @setEvalBranchQuota(20000);
    var table: FaceEdgeAdjTable = undefined;
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        var e: usize = 0;
        while (e < 3) : (e += 1) {
            const face = icosa_face_vertices[f];
            const a = face[(e + 1) % 3];
            const b = face[(e + 2) % 3];
            var found = false;
            var nf: usize = 0;
            while (nf < 20) : (nf += 1) {
                if (nf == f) continue;
                const nface = icosa_face_vertices[nf];
                var ne: usize = 0;
                while (ne < 3) : (ne += 1) {
                    const na = nface[(ne + 1) % 3];
                    const nb = nface[(ne + 2) % 3];
                    if ((na == a and nb == b) or (na == b and nb == a)) {
                        table[f][e] = .{
                            .nf = nf,
                            .ne = ne,
                            .rev = (na == b and nb == a),
                        };
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (!found) @compileError("Face adjacency missing");
        }
    }
    return table;
}

fn maxNeighbors(is_penta: bool) u8 {
    return if (is_penta) 5 else 6;
}

fn hasNeighbor(bt: *const BuilderTile, candidate: usize) bool {
    var i: usize = 0;
    while (i < bt.neighbor_count) : (i += 1) {
        if (bt.tile.neighbors[i] == candidate) return true;
    }
    return false;
}

fn addNeighborOneWay(tiles: []BuilderTile, idx: usize, neighbor: usize) bool {
    var entry = &tiles[idx];
    if (hasNeighbor(entry, neighbor)) return false;
    const cap = maxNeighbors(entry.tile.is_pentagon);
    if (entry.neighbor_count >= cap) return false;
    entry.tile.neighbors[entry.neighbor_count] = @intCast(neighbor);
    entry.neighbor_count += 1;
    return true;
}

fn addNeighborPair(tiles: []BuilderTile, a: usize, b: usize) void {
    if (a == b) return;
    const added_a = addNeighborOneWay(tiles, a, b);
    const added_b = addNeighborOneWay(tiles, b, a);
    _ = added_a;
    _ = added_b;
}

fn sortEdgeEntries(entries: []EdgeEntry) void {
    var i: usize = 1;
    while (i < entries.len) : (i += 1) {
        var j = i;
        while (j > 0 and entries[j - 1].key > entries[j].key) : (j -= 1) {
            const tmp = entries[j - 1];
            entries[j - 1] = entries[j];
            entries[j] = tmp;
        }
    }
}

fn stitchSeamNeighbors(tiles: []BuilderTile, edge_lists: *EdgeEntryGrid, edge_len: usize) void {
    if (edge_len == 0) return;
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        var e: usize = 0;
        while (e < 3) : (e += 1) {
            const adj = face_edge_adj[f][e];
            if (adj.nf < f or (adj.nf == f and adj.ne <= e)) continue;
            const list_a = edge_lists[f][e].items;
            const list_b = edge_lists[adj.nf][adj.ne].items;
            std.debug.assert(list_a.len == edge_len and list_b.len == edge_len);
            var pos: usize = 0;
            while (pos < edge_len) : (pos += 1) {
                const a_idx = list_a[pos].tile_index;
                const b_pos = if (adj.rev) (edge_len - 1 - pos) else pos;
                const b_idx = list_b[b_pos].tile_index;
                addNeighborPair(tiles, a_idx, b_idx);
            }
        }
    }
}

fn stitchPentagonFans(tiles: []BuilderTile, corners: *CornerLists) void {
    var gv: usize = 0;
    while (gv < corners.len) : (gv += 1) {
        const list = corners[gv].items;
        if (list.len == 0) continue;
        const entries = list;
        const count = entries.len;
        if (count <= 1) continue;
        var k: usize = 0;
        while (k < count) : (k += 1) {
            const curr = entries[k].tile_index;
            const next = entries[(k + 1) % count].tile_index;
            addNeighborPair(tiles, curr, next);
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

fn dot(a: [3]f64, b: [3]f64) f64 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

fn makeBasis(center: [3]f32) struct { u: [3]f64, v: [3]f64 } {
    const c = normalize(.{ @as(f64, center[0]), @as(f64, center[1]), @as(f64, center[2]) });
    var up: [3]f64 = .{ 0.0, 0.0, 1.0 };
    var u = cross(c, up);
    if (std.math.sqrt(dot(u, u)) < 1e-6) {
        up = .{ 0.0, 1.0, 0.0 };
        u = cross(c, up);
    }
    u = normalize(u);
    const v = normalize(cross(c, u));
    return .{ .u = u, .v = v };
}

fn sortNeighborsClockwise(tiles: []BuilderTile, idx: usize) void {
    var entry = &tiles[idx];
    if (entry.neighbor_count <= 1) {
        var i: usize = entry.neighbor_count;
        while (i < 6) : (i += 1) entry.tile.neighbors[i] = UNSET_NEIGHBOR;
        return;
    }
    const basis = makeBasis(entry.tile.pos_sphere);
    var ids: [6]TileId = undefined;
    var angles: [6]f64 = undefined;
    const count = entry.neighbor_count;
    var i: usize = 0;
    while (i < count) : (i += 1) {
        const nb = entry.tile.neighbors[i];
        ids[i] = nb;
        const neighbor_pos = tiles[nb].tile.pos_sphere;
        const rel = .{
            @as(f64, neighbor_pos[0]) - @as(f64, entry.tile.pos_sphere[0]),
            @as(f64, neighbor_pos[1]) - @as(f64, entry.tile.pos_sphere[1]),
            @as(f64, neighbor_pos[2]) - @as(f64, entry.tile.pos_sphere[2]),
        };
        const x = dot(rel, basis.u);
        const y = dot(rel, basis.v);
        angles[i] = std.math.atan2(y, x);
    }
    var m: usize = 1;
    while (m < count) : (m += 1) {
        var k = m;
        while (k > 0 and angles[k - 1] > angles[k]) : (k -= 1) {
            const tmp_a = angles[k - 1];
            angles[k - 1] = angles[k];
            angles[k] = tmp_a;
            const tmp_id = ids[k - 1];
            ids[k - 1] = ids[k];
            ids[k] = tmp_id;
        }
    }
    i = 0;
    while (i < count) : (i += 1) entry.tile.neighbors[i] = ids[i];
    while (i < 6) : (i += 1) entry.tile.neighbors[i] = UNSET_NEIGHBOR;
}

pub fn build_prebaked(grid: *Grid, alloc: std.mem.Allocator) !BuiltPlanet {
    // Analytic per-face enumeration; no BFS or canonicalization.
    const N: isize = @intCast(grid.size);

    var tiles_list = std.ArrayListUnmanaged(BuilderTile){};
    defer tiles_list.deinit(alloc);

    var map = CoordMap.init(alloc);
    defer map.deinit();

    var corner_lists: CornerLists = undefined;
    {
        var i: usize = 0;
        while (i < corner_lists.len) : (i += 1) corner_lists[i] = .{};
    }
    defer {
        var i: usize = 0;
        while (i < corner_lists.len) : (i += 1) corner_lists[i].deinit(alloc);
    }

    const edge_len: usize = if (grid.size > 0) grid.size - 1 else 0;
    var edge_lists: EdgeEntryGrid = undefined;
    {
        var f_init: usize = 0;
        while (f_init < 20) : (f_init += 1) {
            var e_init: usize = 0;
            while (e_init < 3) : (e_init += 1) {
                edge_lists[f_init][e_init] = .{};
            }
        }
    }
    defer {
        var f_deinit: usize = 0;
        while (f_deinit < 20) : (f_deinit += 1) {
            var e_deinit: usize = 0;
            while (e_deinit < 3) : (e_deinit += 1) {
                edge_lists[f_deinit][e_deinit].deinit(alloc);
            }
        }
    }

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

                const tile_index: usize = tiles_list.items.len;
                try tiles_list.append(alloc, .{
                    .tile = .{
                        .pos_sphere = .{ @floatCast(pos.x), @floatCast(pos.y), @floatCast(pos.z) },
                        .face = @intCast(face),
                        .axial = .{ .q = @intCast(q), .r = @intCast(r) },
                        .neighbors = .{UNSET_NEIGHBOR} ** 6,
                        .is_pentagon = is_penta,
                    },
                    .neighbor_count = 0,
                });
                try map.put(makeKey(@intCast(face), @intCast(q), @intCast(r)), @intCast(tile_index));

                if (is_penta) {
                    const face_tri = icosa_face_vertices[face];
                    var local: usize = 0;
                    inline for (0..3) |comp| {
                        if (uvw[comp] == 2 * N) {
                            local = comp;
                            break;
                        }
                    }
                    const global_vertex = face_tri[local] - 1;
                    try corner_lists[global_vertex].append(alloc, .{
                        .face = face,
                        .tile_index = tile_index,
                        .vertex = @intCast(local),
                    });
                } else if (edge_len > 0) {
                    var edge_opt: ?usize = null;
                    inline for (0..3) |comp| {
                        if (uvw[comp] == 0) {
                            edge_opt = comp;
                            break;
                        }
                    }
                    if (edge_opt) |edge_idx| {
                        const key_component = uvw[(edge_idx + 1) % 3];
                        try edge_lists[face][edge_idx].append(alloc, .{
                            .tile_index = tile_index,
                            .key = key_component,
                        });
                    }
                }
            }
        }
    }

    const builder_tiles = try tiles_list.toOwnedSlice(alloc);
    defer alloc.free(builder_tiles);

    // Sanity: count per-face tiles and pentagons.
    var per_face = [_]usize{0} ** 20;
    var pentas: usize = 0;
    for (builder_tiles) |bt| {
        per_face[bt.tile.face] += 1;
        if (bt.tile.is_pentagon) pentas += 1;
    }
    std.debug.print("PREBAKE: per-face counts {any}\n", .{per_face});
    std.debug.print("PREBAKE: pentagons {d}\n", .{pentas});

    // Fill neighbor references.
    var ti: usize = 0;
    while (ti < builder_tiles.len) : (ti += 1) {
        const t = builder_tiles[ti].tile;
        for (DIRS) |dir| {
            const nq = t.axial.q + dir[0];
            const nr = t.axial.r + dir[1];
            const key = makeKey(t.face, nq, nr);
            if (map.get(key)) |nb_id| {
                addNeighborPair(builder_tiles, ti, nb_id);
            }
        }
    }

    if (edge_len > 0) {
        var f_sort: usize = 0;
        while (f_sort < 20) : (f_sort += 1) {
            var e_sort: usize = 0;
            while (e_sort < 3) : (e_sort += 1) {
                const list = &edge_lists[f_sort][e_sort];
                if (list.items.len != edge_len) {
                    std.debug.panic(
                        "edge mismatch face={d} edge={d} len={d} expected={d}\n",
                        .{ f_sort, e_sort, list.items.len, edge_len },
                    );
                }
                sortEdgeEntries(list.items);
            }
        }
        stitchSeamNeighbors(builder_tiles, &edge_lists, edge_len);
    }

    stitchPentagonFans(builder_tiles, &corner_lists);

    var sort_idx: usize = 0;
    while (sort_idx < builder_tiles.len) : (sort_idx += 1) {
        sortNeighborsClockwise(builder_tiles, sort_idx);
    }

    var tiles = try alloc.alloc(Tile, builder_tiles.len);
    errdefer alloc.free(tiles);
    for (builder_tiles, 0..) |bt, idx| tiles[idx] = bt.tile;

    const mesh = try buildMeshBuffers(alloc, tiles);
    return .{
        .tiles = tiles,
        .vertices = mesh.vertices,
        .elevations = mesh.elevations,
        .indices = mesh.indices,
        .edges = mesh.edges,
    };
}

const MeshBuffers = struct {
    vertices: []f32,
    elevations: []f32,
    indices: []u32,
    edges: []u32,
};

fn quantizeKey(pos: [3]f32) u128 {
    const scale: f64 = 1_000_000.0;
    var key: u128 = 0;
    inline for (0..3) |i| {
        const value = std.math.round(@as(f64, pos[i]) * scale);
        const quant: i32 = @intFromFloat(value);
        const quant_bits: u32 = @bitCast(quant);
        const bits: u64 = @intCast(quant_bits);
        const shift: u7 = @intCast((2 - i) * 32);
        key |= @as(u128, bits) << shift;
    }
    return key;
}

fn buildMeshBuffers(alloc: std.mem.Allocator, tiles: []const Tile) !MeshBuffers {
    var vertex_map = std.AutoHashMap(u128, u32).init(alloc);
    defer vertex_map.deinit();

    var vertex_data = std.ArrayListUnmanaged(f32){};
    defer vertex_data.deinit(alloc);

    var tile_to_vertex = try alloc.alloc(u32, tiles.len);
    defer alloc.free(tile_to_vertex);

    for (tiles, 0..) |tile, idx| {
        const key = quantizeKey(tile.pos_sphere);
        if (vertex_map.get(key)) |vid| {
            tile_to_vertex[idx] = vid;
            continue;
        }

        const vid: u32 = @intCast(vertex_data.items.len / 3);
        try vertex_data.appendSlice(alloc, &tile.pos_sphere);
        try vertex_map.put(key, vid);
        tile_to_vertex[idx] = vid;
    }

    const vertices = try vertex_data.toOwnedSlice(alloc);
    errdefer alloc.free(vertices);

    const vertex_count: usize = vertices.len / 3;
    const elevations = try alloc.alloc(f32, vertex_count);
    errdefer alloc.free(elevations);
    @memset(elevations, 0);

    var tri_list = std.ArrayListUnmanaged(u32){};
    defer tri_list.deinit(alloc);

    var edge_map = std.AutoHashMap(u64, u8).init(alloc);
    defer edge_map.deinit();
    var edge_list = std.ArrayListUnmanaged(u32){};
    defer edge_list.deinit(alloc);

    var neighbor_buf: [6]usize = undefined;

    var i: usize = 0;
    while (i < tiles.len) : (i += 1) {
        const tile = tiles[i];
        const a_vid = tile_to_vertex[i];

        var neighbor_count: usize = 0;
        for (tile.neighbors) |nb| {
            if (nb == UNSET_NEIGHBOR) continue;
            const nb_idx: usize = @intCast(nb);
            if (nb_idx >= tiles.len) continue;
            var seen = false;
            var s: usize = 0;
            while (s < neighbor_count) : (s += 1) {
                if (neighbor_buf[s] == nb_idx) {
                    seen = true;
                    break;
                }
            }
            if (!seen and neighbor_count < neighbor_buf.len) {
                neighbor_buf[neighbor_count] = nb_idx;
                neighbor_count += 1;
            }
        }

        if (neighbor_count >= 2) {
            var j: usize = 0;
            while (j < neighbor_count) : (j += 1) {
                const nb_b = neighbor_buf[j];
                const nb_c = neighbor_buf[(j + 1) % neighbor_count];
                if (!(i < nb_b and i < nb_c)) continue;

                const b_vid = tile_to_vertex[nb_b];
                const c_vid = tile_to_vertex[nb_c];
                if (a_vid == b_vid or a_vid == c_vid or b_vid == c_vid) continue;

                try tri_list.appendSlice(alloc, &.{ a_vid, b_vid, c_vid });
            }
        }

        var e: usize = 0;
        while (e < neighbor_count) : (e += 1) {
            const nb_idx = neighbor_buf[e];
            const b_vid = tile_to_vertex[nb_idx];
            if (a_vid == b_vid) continue;
            const min = @min(a_vid, b_vid);
            const max = @max(a_vid, b_vid);
            const key = (@as(u64, min) << 32) | @as(u64, max);
            const entry = try edge_map.getOrPut(key);
            if (!entry.found_existing) {
                entry.value_ptr.* = 1;
                try edge_list.appendSlice(alloc, &.{ min, max });
            }
        }
    }

    const indices = try tri_list.toOwnedSlice(alloc);
    errdefer alloc.free(indices);
    const edges = try edge_list.toOwnedSlice(alloc);
    errdefer alloc.free(edges);

    return .{
        .vertices = vertices,
        .elevations = elevations,
        .indices = indices,
        .edges = edges,
    };
}

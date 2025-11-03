const std = @import("std");
const planet = @import("index.zig");
const expect = std.testing.expect;
const expectEqual = std.testing.expectEqual;

// Helper to iterate q,r in [-N..N]
fn eachValid(N: usize, grid: *planet.Grid, f: fn (q: isize, r: isize, g: *planet.Grid) anyerror!void) !void {
    const Ni: isize = @intCast(N);
    var q: isize = -Ni;
    while (q <= Ni) : (q += 1) {
        var r: isize = -Ni;
        while (r <= Ni) : (r += 1) {
            try f(q, r, grid);
        }
    }
}

// Hash roundtrip and face extraction
test "grid: hash roundtrip" {
    var g = try planet.Grid.init(std.testing.allocator, 8);
    defer g.deinit();
    for (0..20) |face| {
        const Ni: isize = @intCast(8);
        var q: isize = -Ni;
        while (q <= Ni) : (q += 1) {
            var r: isize = -Ni;
            while (r <= Ni) : (r += 1) {
                if (!g.isValid(q, r)) continue;
                const h = planet.Grid.hashCoord(q, r, face);
                const t = planet.Grid.unhashCoord(h);
                try expectEqual(q, t.q);
                try expectEqual(r, t.r);
                try expectEqual(@as(usize, face), t.face);
            }
        }
    }
}

// Valid/edge/penta partitioning
test "grid: valid/edge/penta partition" {
    var g = try planet.Grid.init(std.testing.allocator, 6);
    defer g.deinit();
    try eachValid(6, &g, struct {
        fn run(q: isize, r: isize, grid: *planet.Grid) !void {
            const v = grid.isValid(q, r);
            const e = grid.isEdge(q, r);
            const p = grid.isPenta(q, r);
            if (!v) return;
            // Exactly one of inner, edge_nonpenta, penta
            const inner = v and !e and !p;
            const edge_nonpenta = v and e and !p;
            const is_penta = v and p;
            try expect(@intFromBool(inner) + @intFromBool(edge_nonpenta) + @intFromBool(is_penta) == 1);
        }
    }.run);
}

// Edge position sanity and mirror
test "grid: edge pos and mirror" {
    var g = try planet.Grid.init(std.testing.allocator, 7);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();
    const N = g.size;
    if (N <= 1) return; // nothing to test
    for (0..20) |face| {
        const Ni: isize = @intCast(N);
        var q: isize = -Ni;
        while (q <= Ni) : (q += 1) {
            var r: isize = -Ni;
            while (r <= Ni) : (r += 1) {
                if (!g.isValid(q, r) or !g.isEdge(q, r) or g.isPenta(q, r)) continue;
                const ei_opt = g.testGetEdgeIndex(q, r);
                try expect(ei_opt != null);
                const ei = ei_opt.?;
                const pos = g.testGetEdgePosition(q, r, ei);
                try expect(pos < N - 1);
                const pair_opt = g.testEdgeAliasPair(q, r, face);
                try expect(pair_opt != null);
                const pair = pair_opt.?;
                // Ensure both sides are assigned (not UNSET)
                try expect(pair.a != std.math.maxInt(usize));
                try expect(pair.b != std.math.maxInt(usize));
            }
        }
    }
}

// Pentas deduplicate to exactly 12 indices
test "grid: pentas are 12 unique indices" {
    var g = try planet.Grid.init(std.testing.allocator, 6);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();
    var seen = std.AutoHashMap(usize, void).init(std.testing.allocator);
    defer seen.deinit();
    var count: usize = 0;
    for (0..20) |face| {
        const Ni: isize = @intCast(g.size);
        var q: isize = -Ni;
        while (q <= Ni) : (q += 1) {
            var r: isize = -Ni;
            while (r <= Ni) : (r += 1) {
                if (!g.isValid(q, r) or !g.isPenta(q, r)) continue;
                const h = planet.Grid.hashCoord(q, r, face);
                const idx = m.get(h).?;
                if (!seen.contains(idx)) {
                    try seen.put(idx, {});
                    count += 1;
                }
            }
        }
    }
    try expectEqual(@as(usize, 12), count);
}

// Index consistency after populateIndices
test "grid: indices within bounds and unique mapping" {
    var g = try planet.Grid.init(std.testing.allocator, 3);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();
    var it = m.iterator();
    while (it.next()) |e| {
        try expect(e.value_ptr.* < g.tile_count);
    }
}

// Neighbor ring basic invariants and cross-face steps
test "grid: neighbor rings valid and symmetric-ish" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();
    try g.populateNeighbors(&m);
    for (0..g.tile_count) |i| {
        const nbrs = g.neighbors[i];
        const c = g.coords[i];
        const expected: u8 = if (g.isPenta(c.q, c.r)) 5 else 6;
        try expectEqual(expected, nbrs.count);
        var k: usize = 0;
        while (k < nbrs.count) : (k += 1) try expect(nbrs.ids[k] < g.tile_count);
    }
    // symmetry check (weak)
    for (0..g.tile_count) |a| {
        var ia: usize = 0;
        while (ia < g.neighbors[a].count) : (ia += 1) {
            const b = g.neighbors[a].ids[ia];
            if (b == a) continue; // skip padding/self entries
            var found = false;
            var jb: usize = 0;
            while (jb < g.neighbors[b].count) : (jb += 1) {
                if (g.neighbors[b].ids[jb] == a) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                const ca = g.coords[a];
                const cb = g.coords[b];
                std.debug.print("Asymmetry: a={d} face={d} q={d} r={d} -> b={d} face={d} q={d} r={d}\n", .{ a, ca.face, ca.q, ca.r, b, cb.face, cb.q, cb.r });
                std.debug.print("a nbrs (count={d}): ", .{g.neighbors[a].count});
                var ka: usize = 0;
                while (ka < g.neighbors[a].count) : (ka += 1) {
                    const nb = g.neighbors[a].ids[ka];
                    const cnb = g.coords[nb];
                    std.debug.print("{d}(f={d},q={d},r={d}) ", .{ nb, cnb.face, cnb.q, cnb.r });
                }
                std.debug.print("\n", .{});
                std.debug.print("b nbrs (count={d}): ", .{g.neighbors[b].count});
                var kb: usize = 0;
                while (kb < g.neighbors[b].count) : (kb += 1) {
                    const nb2 = g.neighbors[b].ids[kb];
                    const cnb2 = g.coords[nb2];
                    std.debug.print("{d}(f={d},q={d},r={d}) ", .{ nb2, cnb2.face, cnb2.q, cnb2.r });
                }
                std.debug.print("\n", .{});
                try expect(false);
            }
        }
    }
}

// Tile count formula sanity
test "grid: tile_count formula" {
    const sizes = [_]usize{ 2, 3, 4, 5 };
    for (sizes) |N| {
        var g = try planet.Grid.init(std.testing.allocator, N);
        defer g.deinit();
        var m = try g.populateIndices();
        defer m.deinit();
        // Consistency: tile_count equals number of unique coord->index values
        var unique = std.AutoHashMap(usize, void).init(std.testing.allocator);
        defer unique.deinit();
        var it = m.iterator();
        while (it.next()) |e| {
            _ = try unique.put(e.value_ptr.*, {});
        }
        try expectEqual(unique.count(), g.tile_count);
    }
}

// GenerateMesh placeholder
test "grid: generateMesh placeholder" {
    var g = try planet.Grid.init(std.testing.allocator, 5);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();
    try g.populateNeighbors(&m);
    try g.generateMesh();
}

// Edge alias bijection across all faces/edges and sizes
test "grid: edge alias bijection" {
    const sizes_to_test = [_]usize{ 2, 3, 4, 5, 8 };
    for (sizes_to_test) |N| {
        if (N <= 1) continue;
        const edge_len: usize = N - 1;

        var g = try planet.Grid.init(std.testing.allocator, N);
        defer g.deinit();

        // We don't need to build indices here; this test checks the topology mapping only
        for (0..20) |face_a| {
            for (0..3) |edge_a| {
                const map_fwd = planet.Grid.testNeighborEdgeMap(face_a, edge_a);
                const face_b = map_fwd.nf;
                const edge_b = map_fwd.ne;

                const map_rev = planet.Grid.testNeighborEdgeMap(face_b, edge_b);

                var pos_a: usize = 0;
                while (pos_a < edge_len) : (pos_a += 1) {
                    const pos_b = if (map_fwd.rev) (edge_len - 1) - pos_a else pos_a;
                    const pos_c = if (map_rev.rev) (edge_len - 1) - pos_b else pos_b;
                    try expectEqual(pos_a, pos_c);
                }
            }
        }
    }
}

// Degree (valence) check across sizes
test "grid: degrees are 5 at pentas, 6 elsewhere" {
    const sizes = [_]usize{ 2, 3, 4, 5, 7, 8 };
    for (sizes) |N| {
        var g = try planet.Grid.init(std.testing.allocator, N);
        defer g.deinit();

        var map = try g.populateIndices();
        defer map.deinit();
        try g.populateNeighbors(&map);

        var i: usize = 0;
        while (i < g.tile_count) : (i += 1) {
            const c = g.coords[i];
            const expect_deg: u8 = if (g.isPenta(c.q, c.r)) 5 else 6;
            try expectEqual(expect_deg, g.neighbors[i].count);
        }
    }
}

// Pentagon loop-closure and reciprocity across seams
test "grid: pentagon ring closes consistently across seams" {
    const sizes = [_]usize{ 3, 4, 5, 8 };
    for (sizes) |N| {
        var g = try planet.Grid.init(std.testing.allocator, N);
        defer g.deinit();

        var map = try g.populateIndices();
        defer map.deinit();
        try g.populateNeighbors(&map);

        var i: usize = 0;
        while (i < g.tile_count) : (i += 1) {
            const c = g.coords[i];
            if (!g.isPenta(c.q, c.r)) continue;

            const ring = g.neighbors[i];
            try expectEqual(@as(u8, 5), ring.count);

            // Uniqueness and non-self
            var k: usize = 0;
            while (k < 5) : (k += 1) {
                const n = ring.ids[k];
                try expect(n != i);
                var dup = false;
                var j: usize = 0;
                while (j < k) : (j += 1) {
                    if (ring.ids[j] == n) {
                        dup = true;
                    }
                }
                try expect(!dup);
            }

            // Reciprocity
            k = 0;
            while (k < 5) : (k += 1) {
                const n = ring.ids[k];
                var has_back = false;
                var j: usize = 0;
                while (j < g.neighbors[n].count) : (j += 1) {
                    if (g.neighbors[n].ids[j] == i) {
                        has_back = true;
                        break;
                    }
                }
                try expect(has_back);
            }
        }
    }
}

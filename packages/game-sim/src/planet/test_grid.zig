const std = @import("std");
const planet = @import("index.zig");
const gridmod = @import("grid.zig");
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
                const ei_face_opt = g.testFaceLocalEdgeIndex(face, q, r);
                if (ei_face_opt == null) continue; // not a boundary of this face
                const ei = ei_face_opt.?;
                const pos = g.testGetEdgePositionFace(q, r, face, ei);
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

test "grid: seam assumptions" {
    var g = try planet.Grid.init(std.testing.allocator, 6);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();

    // 1) face_neighbors symmetry and vertex pair equality
    const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        var ei: usize = 0;
        while (ei < 3) : (ei += 1) {
            const map = g.testNeighborEdgeMap(f, ei);
            const nf = map.nf;
            const ne = map.ne;
            const A = gridmod.icosa_face_vertices[f];
            const B = gridmod.icosa_face_vertices[nf];
            const a0 = A[pairs[ei][0]];
            const a1 = A[pairs[ei][1]];
            const b0 = B[pairs[ne][0]];
            const b1 = B[pairs[ne][1]];
            try expect((a0 == b0 and a1 == b1) or (a0 == b1 and a1 == b0));
        }
    }

    // 5/6) edge parameterization direction and seam mirroring
    const N: isize = 6;
    f = 0;
    while (f < 20) : (f += 1) {
        var ei: usize = 0;
        while (ei < 3) : (ei += 1) {
            // probe interior positions along edge
            var pos: usize = 0;
            while (pos < @as(usize, @intCast(N - 1))) : (pos += 1) {
                if (pos >= @as(usize, @intCast(N - 1))) break;
                // build uvw on edge (face-local): c=0, b = N+1+2*pos, a = 3N - b
                const a_local = pairs[ei][0];
                const b_local = pairs[ei][1];
                const c_local: usize = 3 - a_local - b_local;
                var uvw_f: [3]isize = .{ 0, 0, 0 };
                uvw_f[c_local] = 0;
                uvw_f[b_local] = (N + 1) + 2 * @as(isize, @intCast(pos));
                uvw_f[a_local] = 3 * N - uvw_f[b_local];
                const qr = g.testFromBaryFace(f, uvw_f);
                if (!g.isEdge(qr.q, qr.r) or g.isPenta(qr.q, qr.r)) continue;
                const pos_check = g.testGetEdgePositionFace(qr.q, qr.r, f, ei);
                try expectEqual(@as(usize, pos), pos_check);
                const uvw_back = g.testToBaryFace(f, qr.q, qr.r);
                try expect((uvw_back[b_local] & 1) == 1);
                try expect(uvw_back[b_local] >= (N + 1) and uvw_back[b_local] <= (2 * N - 1));

                const map = g.testNeighborEdgeMap(f, ei);
                // Stepping across: try all dirs and pick any that changes face
                var dir: usize = 0;
                var checked = false;
                while (dir < 6) : (dir += 1) {
                    const d = g.testStepAcrossOrIn(qr.q, qr.r, f, @intCast(dir));
                    if (d.face != f) {
                        const pos_nf = g.testGetEdgePositionFace(d.q, d.r, map.nf, map.ne);
                        const L: usize = @intCast(N - 1);
                        const rev_pos: usize = if (map.rev) (L - 1) - pos else pos;
                        try expectEqual(rev_pos, pos_nf);
                        checked = true;
                        break;
                    }
                }
                try expect(checked);
            }
        }
    }
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
                const map_fwd = g.testNeighborEdgeMap(face_a, edge_a);
                const face_b = map_fwd.nf;
                const edge_b = map_fwd.ne;

                const map_rev = g.testNeighborEdgeMap(face_b, edge_b);
                // Roundtrip should return original (face, edge)
                try expectEqual(@as(usize, face_a), map_rev.nf);
                try expectEqual(@as(usize, edge_a), map_rev.ne);

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

// Edge alias equality: both sides of an aliased edge map to the same index
test "grid: edge alias indices equal across seam" {
    const sizes = [_]usize{ 3, 4, 5, 8 };
    for (sizes) |N| {
        if (N <= 1) continue;
        var g = try planet.Grid.init(std.testing.allocator, N);
        defer g.deinit();
        var m = try g.populateIndices();
        defer m.deinit();

        for (0..20) |face| {
            const Ni: isize = @intCast(N);
            var q: isize = -Ni;
            while (q <= Ni) : (q += 1) {
                var r: isize = -Ni;
                while (r <= Ni) : (r += 1) {
                    if (!g.isValid(q, r) or !g.isEdge(q, r) or g.isPenta(q, r)) continue;
                    const ei_opt = g.testGetEdgeIndex(q, r);
                    if (ei_opt == null) continue;
                    const pair_opt = g.testEdgeAliasPair(q, r, face);
                    try expect(pair_opt != null);
                    const pair = pair_opt.?;
                    try expect(pair.a == pair.b);
                }
            }
        }
    }
}

// Invariant: every stepAcrossOrIn output is indexable
test "grid: step outputs are always indexable" {
    var g = try planet.Grid.init(std.testing.allocator, 6);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();

    for (0..20) |face| {
        const Ni: isize = @intCast(g.size);
        var q: isize = -Ni;
        while (q <= Ni) : (q += 1) {
            var r: isize = -Ni;
            while (r <= Ni) : (r += 1) {
                if (!g.isValid(q, r)) continue;
                var d: u8 = 0;
                while (d < 6) : (d += 1) {
                    const out = g.testStepAcrossOrIn(q, r, face, d);
                    try expect(g.isValid(out.q, out.r));
                    const h = planet.Grid.hashCoord(out.q, out.r, out.face);
                    if (m.get(h) == null) {
                        const uvw_in = g.testToBaryFace(face, q, r);
                        const uvw_out = g.testToBaryFace(out.face, out.q, out.r);
                        std.debug.print("MISS at face={d} q={d} r={d} dir={d} -> nf={d} oq={d} or={d}\n", .{ face, q, r, d, out.face, out.q, out.r });
                        std.debug.print("  uvw_in =({d},{d},{d})  uvw_out=({d},{d},{d})\n", .{ uvw_in[0], uvw_in[1], uvw_in[2], uvw_out[0], uvw_out[1], uvw_out[2] });
                        std.debug.print("  axes[f]=[{d},{d},{d}] axes[nf]=[{d},{d},{d}]\n", .{
                            g.face_axes[face][0],     g.face_axes[face][1],     g.face_axes[face][2],
                            g.face_axes[out.face][0], g.face_axes[out.face][1], g.face_axes[out.face][2],
                        });
                        return error.TestUnexpectedResult;
                    }
                }
            }
        }
    }
}

// Neighbor uniqueness and degree invariants across tile classes
test "grid: neighbor uniqueness and degree invariants" {
    const sizes = [_]usize{ 3, 4, 6, 8 };
    for (sizes) |N| {
        var g = try planet.Grid.init(std.testing.allocator, N);
        defer g.deinit();
        var m = try g.populateIndices();
        defer m.deinit();
        try g.populateNeighbors(&m);

        var i: usize = 0;
        while (i < g.tile_count) : (i += 1) {
            const c = g.coords[i];
            const is_edge = g.isEdge(c.q, c.r);
            const is_penta = g.isPenta(c.q, c.r);
            const expected: u8 = if (is_penta) 5 else 6;
            const ring = g.neighbors[i];
            try expectEqual(expected, ring.count);
            // Uniqueness and non-self
            var a: usize = 0;
            while (a < ring.count) : (a += 1) {
                const va = ring.ids[a];
                try expect(va != i);
                var b: usize = a + 1;
                while (b < ring.count) : (b += 1) {
                    try expect(ring.ids[b] != va);
                }
            }
            // Edge-interior tiles must still have 6 unique neighbors
            if (is_edge and !is_penta) {
                try expectEqual(@as(u8, 6), ring.count);
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

// Projection central-symmetry property: proj(C+d) + proj(C-d) == 2*C
test "grid: enforceHexWindow symmetry" {
    const N: isize = 8;
    const center: [3]isize = .{ N, N, N };

    const displacements = [_][3]isize{
        .{ 1, 0, -1 },
        .{ 0, 1, -1 },
        .{ 1, -1, 0 },
        .{ -1, 0, 1 },
        .{ 0, -1, 1 },
        .{ -1, 1, 0 },
        .{ 2, -1, -1 },
        .{ N, -N, 0 },
        .{ N + 1, -N, -1 },
        .{ N, 1, -(N + 1) },
    };

    for (displacements) |d| {
        const p_fwd: [3]isize = .{ center[0] + d[0], center[1] + d[1], center[2] + d[2] };
        const p_bwd: [3]isize = .{ center[0] - d[0], center[1] - d[1], center[2] - d[2] };

        const proj_fwd = gridmod.Grid.testEnforceHexWindow(p_fwd, N);
        const proj_bwd = gridmod.Grid.testEnforceHexWindow(p_bwd, N);

        const sum: [3]isize = .{ proj_fwd[0] + proj_bwd[0], proj_fwd[1] + proj_bwd[1], proj_fwd[2] + proj_bwd[2] };
        const expected: [3]isize = .{ 2 * N, 2 * N, 2 * N };

        if (!(sum[0] == expected[0] and sum[1] == expected[1] and sum[2] == expected[2])) {
            std.debug.print("Symmetry failed N={d} d=({d},{d},{d}) sum=({d},{d},{d})\n", .{ N, d[0], d[1], d[2], sum[0], sum[1], sum[2] });
            try expect(false);
        }
    }
}

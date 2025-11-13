const std = @import("std");
const planet = @import("index.zig");
const gridmod = @import("grid.zig");
const expect = std.testing.expect;
const expectEqual = std.testing.expectEqual;

const builtin = @import("builtin");

const faces_to_check = [_]usize{ 0, 1, 2, 8, 12 };

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

test "grid: debug asym pair step (f8,q2,r0) <-> (f17,q-1,r2)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const A = gridmod.Coord{ .q = 2, .r = 0, .face = 8 };
    const B = gridmod.Coord{ .q = -1, .r = 2, .face = 17 };
    std.debug.print("FROM A (f{d},q{d},r{d}):\n", .{ A.face, A.q, A.r });
    var d: u8 = 0;
    while (d < 6) : (d += 1) {
        const nb = g.testStep(A.q, A.r, A.face, d);
        std.debug.print("  dir={d} -> (f{d},q{d},r{d})\n", .{ d, nb.face, nb.q, nb.r });
    }
    std.debug.print("FROM B (f{d},q{d},r{d}):\n", .{ B.face, B.q, B.r });
    d = 0;
    while (d < 6) : (d += 1) {
        const nb = g.testStep(B.q, B.r, B.face, d);
        std.debug.print("  dir={d} -> (f{d},q{d},r{d})\n", .{ d, nb.face, nb.q, nb.r });
    }
}

test "grid: debug asym A-side stepAcrossOrIn (f8,q2,r0) all dirs" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const A = gridmod.Coord{ .q = 2, .r = 0, .face = 8 };
    std.debug.print("A-SIDE stepAcrossOrIn FROM (f{d},q{d},r{d}):\n", .{ A.face, A.q, A.r });
    var d: u8 = 0;
    while (d < 6) : (d += 1) {
        const nb = g.testStepAcrossOrIn(A.q, A.r, A.face, d);
        std.debug.print("  dir={d} -> (f{d},q{d},r{d})\n", .{ d, nb.face, nb.q, nb.r });
    }
}

test "grid: bary inverse oracle vs toBaryFace for A across neighbors" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const N: isize = @intCast(g.size);
    const Aface: usize = 8;
    const Aq: isize = 2;
    const Ar: isize = 0;
    const posA = g.faceCenterPosition(Aq, Ar, Aface);
    // toBaryFace (model) on neighbor faces
    var e: usize = 0;
    while (e < 3) : (e += 1) {
        const nf = g.face_neighbors[Aface][e][0] - 1;
        const model_uvw = g.testToBaryFace(nf, Aq, Ar);
        const gt = gridmod.Grid.baryFromPointOnFace(posA.x, posA.y, posA.z, nf);
        // scale ground-truth to the integer scheme (sum=3N)
        const scale: f64 = 3.0 * @as(f64, @floatFromInt(N));
        const gt_scaled = .{ gt[0] * scale, gt[1] * scale, gt[2] * scale };
        std.debug.print(
            "BARY_ORACLE: A=(f{d},q{d},r{d}) -> nf={d} | model=({d},{d},{d}) gt_scaled=({d:.3},{d:.3},{d:.3})\n",
            .{ Aface, Aq, Ar, nf, model_uvw[0], model_uvw[1], model_uvw[2], gt_scaled[0], gt_scaled[1], gt_scaled[2] },
        );
    }
}

fn eachValidSampled(
    N: usize,
    grid: *planet.Grid,
    max_samples: usize,
    f: fn (q: isize, r: isize, g: *planet.Grid) anyerror!void,
) !void {
    const Ni: isize = @intCast(N);
    var taken: usize = 0;
    // thin cross
    var q: isize = -Ni;
    while (q <= Ni) : (q += 1) {
        const r: isize = 0;
        if (grid.isValid(q, r)) {
            try f(q, r, grid);
            taken += 1;
            if (taken >= max_samples) return;
        }
    }
    var r2: isize = -Ni;
    while (r2 <= Ni) : (r2 += 1) {
        const q2: isize = 0;
        if (grid.isValid(q2, r2)) {
            try f(q2, r2, grid);
            taken += 1;
            if (taken >= max_samples) return;
        }
    }
    // sprinkle interior
    var step: isize = if (Ni <= 2) 1 else (@divTrunc(Ni, 2));
    if (step == 0) step = 1;
    q = -Ni;
    while (q <= Ni and taken < max_samples) : (q += step) {
        var r: isize = -Ni;
        while (r <= Ni and taken < max_samples) : (r += step) {
            if (!grid.isValid(q, r)) continue;
            try f(q, r, grid);
            taken += 1;
        }
    }
}

// Hash roundtrip and face extraction
test "grid: hash roundtrip" {
    var g = try planet.Grid.init(std.testing.allocator, 3);
    defer g.deinit();
    for (faces_to_check) |face| {
        const Ni: isize = @intCast(3);
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
    var g = try planet.Grid.init(std.testing.allocator, 3);
    defer g.deinit();
    try eachValidSampled(3, &g, 200, struct {
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
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();
    const N = g.size;
    if (N <= 1) return; // nothing to test
    for (faces_to_check) |face| {
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
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();
    // Count unique pentagon indices via global canonicalization (gv->idx)
    var seen = std.AutoHashMap(usize, void).init(std.testing.allocator);
    defer seen.deinit();
    var count: usize = 0;
    inline for (0..12) |gv| {
        const idx = g.penta_indices[gv];
        if (idx != 0 and !seen.contains(idx)) {
            try seen.put(idx, {});
            count += 1;
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
    var checked: usize = 0;
    const max_checked: usize = 500;
    var i: usize = 0;
    while (i < g.tile_count and checked < max_checked) : (i += 1) {
        const nbrs = g.neighbors[i];
        const c = g.coords[i];
        const expected: u8 = if (g.isPentaAt(c.face, c.q, c.r)) 5 else 6;
        try expectEqual(expected, nbrs.count);
        var k: usize = 0;
        while (k < nbrs.count) : (k += 1) try expect(nbrs.ids[k] < g.tile_count);
        checked += 1;
    }
    // symmetry (sampled)
    checked = 0;
    var a: usize = 0;
    while (a < g.tile_count and checked < max_checked) : (a += 1) {
        var ia: usize = 0;
        while (ia < g.neighbors[a].count) : (ia += 1) {
            const b = g.neighbors[a].ids[ia];
            if (b == a) continue;
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
                std.debug.print("ASYM: a(i{d}) (f{d},q{d},r{d}) -> b(i{d}) (f{d},q{d},r{d}) not reciprocal\n", .{
                    a, ca.face, ca.q, ca.r, b, cb.face, cb.q, cb.r,
                });
            }
            try expect(found);
        }
        checked += 1;
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

test "seam tables consistent" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    try g.assertSeamTablesConsistent();
}

test "finish invertibility holds on all faces" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const N: isize = @intCast(g.size);
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        var q: isize = -N;
        while (q <= N) : (q += 1) {
            var r: isize = -N;
            while (r <= N) : (r += 1) {
                const c1 = g.canonicalCoord(f, q, r);
                const c2 = g.canonicalCoord(c1.face, c1.q, c1.r);
                try expectEqual(c1.face, c2.face);
                try expectEqual(c1.q, c2.q);
                try expectEqual(c1.r, c2.r);
                // Interior points should be no-ops
                const uvw_face = g.testToBaryFace(f, q, r);
                const interior = (uvw_face[0] > 0 and uvw_face[1] > 0 and uvw_face[2] > 0 and
                    uvw_face[0] < N and uvw_face[1] < N and uvw_face[2] < N);
                if (interior) {
                    try expectEqual(f, c1.face);
                    try expectEqual(q, c1.q);
                    try expectEqual(r, c1.r);
                }
            }
        }
    }
}

test "edge owner canonicalization" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const N: isize = @intCast(g.size);
    const R: isize = 3 * N;
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        var e: usize = 0;
        while (e < 3) : (e += 1) {
            const nf = g.faceNeighborFace(f, e);
            const owner: usize = if (f < nf) f else nf;
            // sample two positions along the edge excluding corners
            const samples = [_]isize{ 1, R - 1 };
            for (samples) |t| {
                var lab: [3]isize = .{ 0, 0, 0 };
                lab[(e + 1) % 3] = t;
                lab[(e + 2) % 3] = R - t;
                const uvw_act = g.testLabelToActual(f, lab);
                const qr = g.testFromBaryFace(f, uvw_act);
                const c = g.canonicalCoord(f, qr.q, qr.r);
                try expectEqual(owner, c.face);
            }
        }
    }
}

test "vertex owner canonicalization" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const N: isize = @intCast(g.size);
    const R: isize = 3 * N;
    // For each global vertex label, craft lab with two zeros and assert owner
    var v: u8 = 0;
    while (v < 12) : (v += 1) {
        // Find a face that contains label v
        var found_face: usize = 20;
        var local_k: usize = 3;
        var f: usize = 0;
        while (f < 20 and found_face == 20) : (f += 1) {
            var k: usize = 0;
            while (k < 3) : (k += 1) {
                if (g.testGetFaceLabel(f, k) == v) {
                    found_face = f;
                    local_k = k;
                    break;
                }
            }
        }
        try expect(found_face < 20);
        var lab: [3]isize = .{ 0, 0, 0 };
        lab[local_k] = R;
        const uvw_act = g.testLabelToActual(found_face, lab);
        const qr = g.testFromBaryFace(found_face, uvw_act);
        const c = g.canonicalCoord(found_face, qr.q, qr.r);
        try expectEqual(@as(usize, g.testVertexOwner(v)), c.face);
    }
}

test "no seam hop on boundary zero" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    // Use the failing dump's label-basis coords: (7,5,0) with R=12 on face 0
    const uvw_lbl = [3]isize{ 7, 5, 0 };
    const uvw_act = g.testLabelToActual(0, uvw_lbl);
    const pick = g.testPickViolatedEdge(0, uvw_act);
    try expect(pick == null);
}

test "canonical interior idempotent (random samples per face)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const N: isize = @intCast(g.size);
    const R: isize = 3 * N;
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        var a_i: isize = 1;
        var count: usize = 0;
        const max_a: isize = if (R - 2 > 20) 20 else (R - 2);
        while (a_i <= max_a and count < 20) : (a_i += 1) {
            const a: isize = a_i;
            const b: isize = 1; // ensures c = R - a - 1 >= 1 for a <= R-2
            const c: isize = R - a - b;
            const lab = [3]isize{ a, b, c };
            const uvw_act = g.testLabelToActual(f, lab);
            const qr = g.testFromBaryFace(f, uvw_act);
            const c0 = g.canonicalCoord(f, qr.q, qr.r);
            try expectEqual(f, c0.face);
            try expectEqual(qr.q, c0.q);
            try expectEqual(qr.r, c0.r);
            count += 1;
        }
    }
}

test "edge owner canonicalization (both sides, perm-only boundary)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const N: isize = @intCast(g.size);
    const R: isize = 3 * N;
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        var e: usize = 0;
        while (e < 3) : (e += 1) {
            const nf = g.faceNeighborFace(f, e);
            const owner: usize = if (f < nf) f else nf;
            const ts = [_]isize{ 1, R - 1 };
            for (ts) |t| {
                // side A
                var labA: [3]isize = .{ 0, 0, 0 };
                labA[(e + 1) % 3] = t;
                labA[(e + 2) % 3] = R - t;
                const uvwA = g.testLabelToActual(f, labA);
                const qrA = g.testFromBaryFace(f, uvwA);
                const cA = g.canonicalCoord(f, qrA.q, qrA.r);
                try expectEqual(owner, cA.face);
                // side B (neighbor’s perspective)
                // Build neighbor lab with zero at its edge index; we can just mirror using same pattern with ne
                var labB: [3]isize = .{ 0, 0, 0 };
                // derive ne via neighborEdgeMap helper
                const map = g.testNeighborEdgeMap(f, e);
                const ne_idx = map.ne;
                labB[(ne_idx + 1) % 3] = t;
                labB[(ne_idx + 2) % 3] = R - t;
                const uvwB = g.testLabelToActual(nf, labB);
                const qrB = g.testFromBaryFace(nf, uvwB);
                const cB = g.canonicalCoord(nf, qrB.q, qrB.r);
                try expectEqual(owner, cB.face);
            }
        }
    }
}

test "round-trip invariant via exact reverse delta on interiors" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    // choose interior grid on each face (q,r in {1,2})
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        var q: isize = 1;
        while (q <= 2) : (q += 1) {
            var r: isize = 1;
            while (r <= 2) : (r += 1) {
                var d: u8 = 0;
                while (d < 6) : (d += 1) {
                    const a0 = g.canonicalCoord(f, q, r);
                    const step = g.testExactStep(a0.face, a0.q, a0.r, d);
                    const b = g.canonicalCoord(step.face, step.q, step.r);
                    const du = g.testDirCubeStep(a0.face, d);
                    const rev = [3]isize{ -du[0], -du[1], -du[2] };
                    const d_rev = g.testFindDirByDelta(b.face, rev) orelse return error.Unexpected;
                    const back = g.testExactStep(b.face, b.q, b.r, d_rev);
                    const back_can = g.canonicalCoord(back.face, back.q, back.r);
                    try expectEqual(a0.face, back_can.face);
                    try expectEqual(a0.q, back_can.q);
                    try expectEqual(a0.r, back_can.r);
                }
            }
        }
    }
}
// test "grid: seam assumptions" {
//     var g = try planet.Grid.initWithMode(std.testing.allocator, 3, .full);
//     defer g.deinit();
//     var m = try g.populateIndices();
//     defer m.deinit();

//     // 1) face_neighbors symmetry and vertex pair equality
//     const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
//     var f: usize = 0;
//     while (f < faces_to_check.len) : (f += 1) {
//         const ff = faces_to_check[f];
//         var ei: usize = 0;
//         while (ei < 3) : (ei += 1) {
//             const map = g.testNeighborEdgeMap(ff, ei);
//             const nf = map.nf;
//             const ne = map.ne;
//             const A = gridmod.icosa_face_vertices[ff];
//             const B = gridmod.icosa_face_vertices[nf];
//             const a0 = A[pairs[ei][0]];
//             const a1 = A[pairs[ei][1]];
//             const b0 = B[pairs[ne][0]];
//             const b1 = B[pairs[ne][1]];
//             try expect((a0 == b0 and a1 == b1) or (a0 == b1 and a1 == b0));
//         }
//     }

//     // 5/6) edge parameterization direction and seam mirroring
//     const N: isize = 3;
//     f = 0;
//     while (f < faces_to_check.len) : (f += 1) {
//         const ff = faces_to_check[f];
//         var ei: usize = 0;
//         while (ei < 3) : (ei += 1) {
//             // probe a few interior positions along edge
//             const L: usize = @intCast(N - 1);
//             const probe = [_]usize{ 0, L / 2, if (L > 1) L - 1 else 0 };
//             var pi: usize = 0;
//             while (pi < probe.len) : (pi += 1) {
//                 const pos = probe[pi];
//                 // build uvw on edge (face-local): c=0, b = N+1+2*pos, a = 3N - b
//                 const a_local = pairs[ei][0];
//                 const b_local = pairs[ei][1];
//                 const c_local: usize = 3 - a_local - b_local;
//                 var uvw_f: [3]isize = .{ 0, 0, 0 };
//                 uvw_f[c_local] = 0;
//                 uvw_f[b_local] = (N + 1) + 2 * @as(isize, @intCast(pos));
//                 uvw_f[a_local] = 3 * N - uvw_f[b_local];
//                 const qr = g.testFromBaryFace(ff, uvw_f);
//                 if (!g.isEdge(qr.q, qr.r) or g.isPenta(qr.q, qr.r)) continue;
//                 const pos_check = g.testGetEdgePositionFace(qr.q, qr.r, ff, ei);
//                 try expectEqual(@as(usize, pos), pos_check);
//                 const uvw_back = g.testToBaryFace(ff, qr.q, qr.r);
//                 try expect((uvw_back[b_local] & 1) == 1);
//                 try expect(uvw_back[b_local] >= (N + 1) and uvw_back[b_local] <= (2 * N - 1));

//                 const map = g.testNeighborEdgeMap(ff, ei);
//                 // Deterministic outward direction based on opposite component
//                 const OUT: [3]u8 = .{ 3, 4, 2 }; // (-q), (-r), (-q,+s)
//                 const dir: u8 = OUT[c_local];
//                 const d = g.testStepAcrossOrIn(qr.q, qr.r, ff, dir);
//                 try expect(d.face != ff);
//                 const pos_nf = g.testGetEdgePositionFace(d.q, d.r, map.nf, map.ne);
//                 const rev_pos: usize = if (map.rev) (L - 1) - pos else pos;
//                 try expectEqual(rev_pos, pos_nf);
//             }
//         }
//     }
// }

// // Edge alias bijection across all faces/edges and sizes
// test "grid: edge alias bijection" {
//     const sizes_to_test = [_]usize{ 2, 3 };
//     for (sizes_to_test) |N| {
//         if (N <= 1) continue;
//         const edge_len: usize = N - 1;

//         var g = try planet.Grid.initWithMode(std.testing.allocator, N, .full);
//         defer g.deinit();

//         // We don't need to build indices here; this test checks the topology mapping only
//         for (faces_to_check) |face_a| {
//             for (0..3) |edge_a| {
//                 const map_fwd = g.testNeighborEdgeMap(face_a, edge_a);
//                 const face_b = map_fwd.nf;
//                 const edge_b = map_fwd.ne;

//                 const map_rev = g.testNeighborEdgeMap(face_b, edge_b);
//                 // Roundtrip should return original (face, edge)
//                 try expectEqual(@as(usize, face_a), map_rev.nf);
//                 try expectEqual(@as(usize, edge_a), map_rev.ne);

//                 var pos_a: usize = 0;
//                 while (pos_a < edge_len) : (pos_a += 1) {
//                     const pos_b = if (map_fwd.rev) (edge_len - 1) - pos_a else pos_a;
//                     const pos_c = if (map_rev.rev) (edge_len - 1) - pos_b else pos_b;
//                     try expectEqual(pos_a, pos_c);
//                 }
//             }
//         }
//     }
// }

// // Edge alias equality: both sides of an aliased edge map to the same index
// test "grid: edge alias indices equal across seam" {
//     const sizes = [_]usize{ 3, 4 };
//     for (sizes) |N| {
//         if (N <= 1) continue;
//         var g = try planet.Grid.initWithMode(std.testing.allocator, N, .full);
//         defer g.deinit();
//         var m = try g.populateIndices();
//         defer m.deinit();

//         for (faces_to_check) |face| {
//             const Ni: isize = @intCast(N);
//             var q: isize = -Ni;
//             while (q <= Ni) : (q += 1) {
//                 var r: isize = -Ni;
//                 while (r <= Ni) : (r += 1) {
//                     if (!g.isValid(q, r) or !g.isEdge(q, r) or g.isPenta(q, r)) continue;
//                     const ei_opt = g.testGetEdgeIndex(q, r);
//                     if (ei_opt == null) continue;
//                     const pair_opt = g.testEdgeAliasPair(q, r, face);
//                     try expect(pair_opt != null);
//                     const pair = pair_opt.?;
//                     try expect(pair.a == pair.b);
//                 }
//             }
//         }
//     }
// }

// // Invariant: every stepAcrossOrIn output is indexable
// test "grid: step outputs are always indexable" {
//     var g = try planet.Grid.init(std.testing.allocator, 3);
//     defer g.deinit();
//     var m = try g.populateIndices();
//     defer m.deinit();

//     for (faces_to_check) |face| {
//         const Ni: isize = @intCast(g.size);
//         var q: isize = -Ni;
//         while (q <= Ni) : (q += 1) {
//             var r: isize = -Ni;
//             while (r <= Ni) : (r += 1) {
//                 if (!g.isValid(q, r)) continue;
//                 var d: u8 = 0;
//                 while (d < 6) : (d += 1) {
//                     const out = g.testStepAcrossOrIn(q, r, face, d);
//                     try expect(g.isValid(out.q, out.r));
//                     const h = planet.Grid.hashCoord(out.q, out.r, out.face);
//                     if (m.get(h) == null) {
//                         return error.TestUnexpectedResult;
//                     }
//                 }
//             }
//         }
//     }
// }

// // Neighbor uniqueness and degree invariants across tile classes
// test "grid: neighbor uniqueness and degree invariants" {
//     const sizes = [_]usize{ 3, 4 };
//     for (sizes) |N| {
//         var g = try planet.Grid.init(std.testing.allocator, N);
//         defer g.deinit();
//         var m = try g.populateIndices();
//         defer m.deinit();
//         try g.populateNeighbors(&m);

//         var i: usize = 0;
//         var checked: usize = 0;
//         const max_checked: usize = 600;
//         while (i < g.tile_count and checked < max_checked) : (i += 1) {
//             const c = g.coords[i];
//             const is_edge = g.isEdge(c.q, c.r);
//             const is_penta = g.isPenta(c.q, c.r);
//             const expected: u8 = if (is_penta) 5 else 6;
//             const ring = g.neighbors[i];
//             try expectEqual(expected, ring.count);
//             // Uniqueness and non-self
//             var a: usize = 0;
//             while (a < ring.count) : (a += 1) {
//                 const va = ring.ids[a];
//                 try expect(va != i);
//                 var b: usize = a + 1;
//                 while (b < ring.count) : (b += 1) {
//                     try expect(ring.ids[b] != va);
//                 }
//             }
//             // Edge-interior tiles must still have 6 unique neighbors
//             if (is_edge and !is_penta) {
//                 try expectEqual(@as(u8, 6), ring.count);
//             }
//             checked += 1;
//         }
//     }
// }

// // Degree (valence) check across sizes
// test "grid: degrees are 5 at pentas, 6 elsewhere" {
//     const sizes = [_]usize{ 2, 3, 4 };
//     for (sizes) |N| {
//         var g = try planet.Grid.init(std.testing.allocator, N);
//         defer g.deinit();

//         var map = try g.populateIndices();
//         defer map.deinit();
//         try g.populateNeighbors(&map);

//         var i: usize = 0;
//         while (i < g.tile_count) : (i += 1) {
//             const c = g.coords[i];
//             const expect_deg: u8 = if (g.isPenta(c.q, c.r)) 5 else 6;
//             try expectEqual(expect_deg, g.neighbors[i].count);
//         }
//     }
// }

// // Pentagon loop-closure and reciprocity across seams
// test "grid: pentagon ring closes consistently across seams" {
//     const sizes = [_]usize{ 3, 4 };
//     for (sizes) |N| {
//         var g = try planet.Grid.init(std.testing.allocator, N);
//         defer g.deinit();

//         var map = try g.populateIndices();
//         defer map.deinit();
//         try g.populateNeighbors(&map);

//         var i: usize = 0;
//         while (i < g.tile_count) : (i += 1) {
//             const c = g.coords[i];
//             if (!g.isPenta(c.q, c.r)) continue;

//             const ring = g.neighbors[i];
//             try expectEqual(@as(u8, 5), ring.count);

//             // Uniqueness and non-self
//             var k: usize = 0;
//             while (k < 5) : (k += 1) {
//                 const n = ring.ids[k];
//                 try expect(n != i);
//                 var dup = false;
//                 var j: usize = 0;
//                 while (j < k) : (j += 1) {
//                     if (ring.ids[j] == n) {
//                         dup = true;
//                     }
//                 }
//                 try expect(!dup);
//             }

//             // Reciprocity
//             k = 0;
//             while (k < 5) : (k += 1) {
//                 const n = ring.ids[k];
//                 var has_back = false;
//                 var j: usize = 0;
//                 while (j < g.neighbors[n].count) : (j += 1) {
//                     if (g.neighbors[n].ids[j] == i) {
//                         has_back = true;
//                         break;
//                     }
//                 }
//                 try expect(has_back);
//             }
//         }
//     }
// }

// // Projection central-symmetry property: proj(C+d) + proj(C-d) == 2*C
// test "grid: enforceHexWindow symmetry" {
//     const N: isize = 8;
//     const center: [3]isize = .{ N, N, N };

//     const displacements = [_][3]isize{
//         .{ 1, 0, -1 },
//         .{ 0, 1, -1 },
//         .{ 1, -1, 0 },
//         .{ -1, 0, 1 },
//         .{ 0, -1, 1 },
//         .{ -1, 1, 0 },
//         .{ 2, -1, -1 },
//         .{ N, -N, 0 },
//         .{ N + 1, -N, -1 },
//         .{ N, 1, -(N + 1) },
//     };

//     for (displacements) |d| {
//         const p_fwd: [3]isize = .{ center[0] + d[0], center[1] + d[1], center[2] + d[2] };
//         const p_bwd: [3]isize = .{ center[0] - d[0], center[1] - d[1], center[2] - d[2] };

//         const proj_fwd = gridmod.Grid.testEnforceHexWindow(p_fwd, N);
//         const proj_bwd = gridmod.Grid.testEnforceHexWindow(p_bwd, N);

//         const sum: [3]isize = .{ proj_fwd[0] + proj_bwd[0], proj_fwd[1] + proj_bwd[1], proj_fwd[2] + proj_bwd[2] };
//         const expected: [3]isize = .{ 2 * N, 2 * N, 2 * N };

//         if (!(sum[0] == expected[0] and sum[1] == expected[1] and sum[2] == expected[2])) {
//             std.debug.print("Symmetry failed N={d} d=({d},{d},{d}) sum=({d},{d},{d})\n", .{ N, d[0], d[1], d[2], sum[0], sum[1], sum[2] });
//             try expect(false);
//         }
//     }
// }

test "grid: seam 0<->8 reciprocity and failing tile diagnostics" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();

    // Find edges linking 0 <-> 8
    var e0: usize = 3;
    var e8: usize = 3;
    var ei: usize = 0;
    while (ei < 3) : (ei += 1) {
        const nf = g.faceNeighborFace(0, ei);
        if (nf == 8) {
            e0 = ei;
            break;
        }
    }
    ei = 0;
    while (ei < 3) : (ei += 1) {
        const nf = g.faceNeighborFace(8, ei);
        if (nf == 0) {
            e8 = ei;
            break;
        }
    }
    try expect(e0 < 3);
    try expect(e8 < 3);

    // Check reciprocity of neighbor edge mapping
    const m0 = g.testNeighborEdgeMap(0, e0);
    const m8 = g.testNeighborEdgeMap(8, e8);
    try expectEqual(@as(usize, 8), m0.nf);
    try expectEqual(@as(usize, 0), m8.nf);
    try expectEqual(@as(usize, e8), m0.ne);
    try expectEqual(@as(usize, e0), m8.ne);

    // Build indices to verify existence
    var m = try g.populateIndices();
    defer m.deinit();

    // Failing center (canonicalize start)
    const start_c = g.canonicalCoord(0, -3, -1);
    const fq: isize = start_c.q;
    const fr: isize = start_c.r;
    const ff: usize = start_c.face;

    // dir=2
    const d2_raw = g.testStepAcrossOrIn(fq, fr, ff, 2);
    const d2 = g.canonicalCoord(d2_raw.face, d2_raw.q, d2_raw.r);
    std.debug.print("DIAG: dir=2 -> (f={d},q={d},r={d})\n", .{ d2.face, d2.q, d2.r });
    const h2 = planet.Grid.hashCoord(d2.q, d2.r, d2.face);
    const idx2_opt = m.get(h2);
    const idx2_val: isize = if (idx2_opt) |ix| @as(isize, @intCast(ix)) else -1;
    std.debug.print("  idx={d}\n", .{idx2_val});

    // dir=3
    const d3_raw = g.testStepAcrossOrIn(fq, fr, ff, 3);
    const d3 = g.canonicalCoord(d3_raw.face, d3_raw.q, d3_raw.r);
    std.debug.print("DIAG: dir=3 -> (f={d},q={d},r={d})\n", .{ d3.face, d3.q, d3.r });
    const h3 = planet.Grid.hashCoord(d3.q, d3.r, d3.face);
    const idx3_opt = m.get(h3);
    const idx3_val: isize = if (idx3_opt) |ix| @as(isize, @intCast(ix)) else -1;
    std.debug.print("  idx={d}\n", .{idx3_val});

    // Reverse step reciprocity checks
    const rev2: u8 = (2 + 3) % 6;
    const back2_raw = g.testStepAcrossOrIn(d2.q, d2.r, d2.face, rev2);
    const back2 = g.canonicalCoord(back2_raw.face, back2_raw.q, back2_raw.r);
    std.debug.print("REV: dir=2 back -> (f={d},q={d},r={d})\n", .{ back2.face, back2.q, back2.r });

    const rev3: u8 = (3 + 3) % 6;
    const back3_raw = g.testStepAcrossOrIn(d3.q, d3.r, d3.face, rev3);
    const back3 = g.canonicalCoord(back3_raw.face, back3_raw.q, back3_raw.r);
    std.debug.print("REV: dir=3 back -> (f={d},q={d},r={d})\n", .{ back3.face, back3.q, back3.r });

    // These may not yet hold, but we log and assert presence in index map if available
    if (idx2_opt != null) {
        try expectEqual(@as(isize, fq), back2.q);
        try expectEqual(@as(isize, fr), back2.r);
        try expectEqual(@as(usize, ff), back2.face);
    }
    if (idx3_opt != null) {
        try expectEqual(@as(isize, fq), back3.q);
        try expectEqual(@as(isize, fr), back3.r);
        try expectEqual(@as(usize, ff), back3.face);
    }
}

test "grid: trace from face 0 (-3,4) across all dirs" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    var dir: u8 = 0;
    while (dir < 6) : (dir += 1) {
        const c = g.testStepAcrossOrIn(-3, 4, 0, dir);
        std.debug.print("TRACE_CALL dir={d} -> (f={d},q={d},r={d})\n", .{ dir, c.face, c.q, c.r });
    }
}

test "grid: seam violation dump at (f0,q-3,r4)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    var d: u8 = 0;
    while (d < 6) : (d += 1) {
        g.testSeamDebug(-3, 4, 0, d);
    }
}

test "grid: axes vs geometry on face 0 -> 12 seam" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();

    // Find the edge on face 0 that maps to face 12
    var e0_to_12: usize = 3;
    var e: usize = 0;
    while (e < 3) : (e += 1) {
        if (g.faceNeighborFace(0, e) == 12) {
            e0_to_12 = e;
            break;
        }
    }
    try expect(e0_to_12 < 3);

    // Gather geometry-derived permutation for that edge
    const perm_geom = g.testBuildPermForEdge(0, e0_to_12);

    // Gather axes-derived expected permutation for that seam
    const Pf0 = g.face_axes[0];
    const P12 = g.face_axes[12];
    const perm_axes = g.testExpectedPermFromAxes(Pf0, P12);

    std.debug.print("AXES_GEOM: e0_to_12={d} perm_geom=({d},{d},{d}) perm_axes=({d},{d},{d}) Pf0=({d},{d},{d}) P12=({d},{d},{d})\n", .{ e0_to_12, perm_geom[0], perm_geom[1], perm_geom[2], perm_axes[0], perm_axes[1], perm_axes[2], Pf0[0], Pf0[1], Pf0[2], P12[0], P12[1], P12[2] });

    try expect(perm_geom[0] == perm_axes[0] and perm_geom[1] == perm_axes[1] and perm_geom[2] == perm_axes[2]);

    // Also verify the stored seam perm matches geometry
    const stored = g.face_perm[0][e0_to_12];
    std.debug.print("AXES_GEOM: stored_perm=({d},{d},{d})\n", .{ stored[0], stored[1], stored[2] });
    try expect(stored[0] == perm_geom[0] and stored[1] == perm_geom[1] and stored[2] == perm_geom[2]);
}

test "grid: face 0 neighbor list" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();

    var seen: [3]usize = .{ 999, 999, 999 };
    var e: usize = 0;
    while (e < 3) : (e += 1) {
        const nf = g.faceNeighborFace(0, e);
        std.debug.print("F0 NEI: e={d} -> nf={d}\n", .{ e, nf });
        seen[e] = nf;
        try expect(nf < 20);
    }
    // assert uniqueness
    try expect(!(seen[0] == seen[1] or seen[0] == seen[2] or seen[1] == seen[2]));
}

test "grid: 5-face neighborhood audit around each vertex of face 0" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();

    var vlocal: usize = 0;
    while (vlocal < 3) : (vlocal += 1) {
        const gv1_based = gridmod.icosa_face_vertices[0][vlocal];
        const gv = gv1_based; // 1-based for comparison against face table entries

        // Collect the five faces that include this global vertex
        var faces5: [5]usize = .{ 999, 999, 999, 999, 999 };
        var count: usize = 0;
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            const tri = gridmod.icosa_face_vertices[f];
            if (tri[0] == gv or tri[1] == gv or tri[2] == gv) {
                if (count < 5) faces5[count] = f;
                count += 1;
            }
        }
        std.debug.print("VERTEX {d} (local {d} on face 0): faces sharing it: ", .{ gv, vlocal });
        for (faces5) |ff| std.debug.print("{d} ", .{ff});
        std.debug.print("\n", .{});
        try expect(count == 5);

        // Degree check: within these five, each face should have exactly two neighbors also in the set
        var ok = true;
        for (faces5) |ff| {
            var deg: usize = 0;
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = g.face_neighbors[ff][e][0] - 1;
                var in_set = false;
                for (faces5) |rf| {
                    if (rf == nf) in_set = true;
                }
                if (in_set) deg += 1;
            }
            if (deg != 2) ok = false;
        }
        try expect(ok);

        // Connectivity check within the 5-face induced subgraph (robust vs. greedy cycle)
        var q = std.ArrayListUnmanaged(usize){};
        defer q.deinit(std.testing.allocator);
        var seen = std.AutoHashMap(usize, void).init(std.testing.allocator);
        defer seen.deinit();
        try q.append(std.testing.allocator, faces5[0]);
        try seen.put(faces5[0], {});
        var qi: usize = 0;
        while (qi < q.items.len) : (qi += 1) {
            const cur = q.items[qi];
            var e2: usize = 0;
            while (e2 < 3) : (e2 += 1) {
                const nf = g.face_neighbors[cur][e2][0] - 1;
                var in_set = false;
                for (faces5) |rf| {
                    if (rf == nf) in_set = true;
                }
                if (!in_set) continue;
                if (!seen.contains(nf)) {
                    try seen.put(nf, {});
                    try q.append(std.testing.allocator, nf);
                }
            }
        }
        try expectEqual(@as(usize, 5), seen.count());

        // We don't enforce a specific cycle ordering; degree+connectivity are sufficient here.
    }
}

test "grid: pinpoint broken seam and winding at face 0 vertex0" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };

    const gv = gridmod.icosa_face_vertices[0][0]; // 1-based global vertex id at local vertex 0 of face 0

    // Collect five faces that include gv
    var faces5: [5]usize = .{ 999, 999, 999, 999, 999 };
    var count: usize = 0;
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        const tri = gridmod.icosa_face_vertices[f];
        if (tri[0] == gv or tri[1] == gv or tri[2] == gv) {
            if (count < 5) faces5[count] = f;
            count += 1;
        }
    }
    std.debug.print("VTX0 gv={d}: faces=", .{gv});
    for (faces5) |ff| std.debug.print("{d} ", .{ff});
    std.debug.print("\n", .{});
    try expect(count == 5);

    // For each face in the set, enumerate edges incident to gv and check neighbor + reciprocity
    for (faces5) |ff| {
        const tri = gridmod.icosa_face_vertices[ff];
        var e: usize = 0;
        while (e < 3) : (e += 1) {
            const a = tri[pairs[e][0]];
            const b = tri[pairs[e][1]];
            const incident = (a == gv or b == gv);
            if (!incident) continue;
            const nf = g.face_neighbors[ff][e][0] - 1;
            const ne = g.face_neighbors[ff][e][1] - 1;
            var in_set = false;
            var si: usize = 0;
            while (si < 5) : (si += 1) {
                if (faces5[si] == nf) {
                    in_set = true;
                    break;
                }
            }
            // Locate the reciprocal edge on nf -> ff
            var back_e: usize = 3;
            var k: usize = 0;
            while (k < 3) : (k += 1) {
                if (g.face_neighbors[nf][k][0] - 1 == ff) {
                    back_e = k;
                    break;
                }
            }
            const tri_nf = gridmod.icosa_face_vertices[nf];
            const a_b = if (back_e < 3) tri_nf[pairs[back_e][0]] else 0;
            const b_b = if (back_e < 3) tri_nf[pairs[back_e][1]] else 0;
            const rec_ok = (back_e < 3) and (g.face_neighbors[nf][back_e][0] - 1 == ff);
            const share_b = (a_b == gv or b_b == gv);
            std.debug.print(
                "SEAM_CHECK: f={d} e={d} (ga,gb)=({d},{d}) -> nf={d} ne={d} in_set={any} rec_ok={any} nf_edge_gv_ok={any}\n",
                .{ ff, e, a, b, nf, ne, in_set, rec_ok, share_b },
            );
        }
    }

    // Winding consistency across the five faces
    const v = gridmod.icosa_vertices;
    const sign_ref: f64 = blk: {
        const f0 = faces5[0];
        const tri0 = gridmod.icosa_face_vertices[f0];
        const ia = tri0[0] - 1;
        const ib = tri0[1] - 1;
        const ic = tri0[2] - 1;
        const a = v[ia];
        const b = v[ib];
        const c = v[ic];
        const ab = .{ .x = b.x - a.x, .y = b.y - a.y, .z = b.z - a.z };
        const ac = .{ .x = c.x - a.x, .y = c.y - a.y, .z = c.z - a.z };
        const n = .{ .x = ab.y * ac.z - ab.z * ac.y, .y = ab.z * ac.x - ab.x * ac.z, .z = ab.x * ac.y - ab.y * ac.x };
        const dot = n.x * a.x + n.y * a.y + n.z * a.z;
        break :blk dot;
    };
    var all_same = true;
    for (faces5) |ff| {
        const tri = gridmod.icosa_face_vertices[ff];
        const ia = tri[0] - 1;
        const ib = tri[1] - 1;
        const ic = tri[2] - 1;
        const a = v[ia];
        const b = v[ib];
        const c = v[ic];
        const ab = .{ .x = b.x - a.x, .y = b.y - a.y, .z = b.z - a.z };
        const ac = .{ .x = c.x - a.x, .y = c.y - a.y, .z = c.z - a.z };
        const n = .{ .x = ab.y * ac.z - ab.z * ac.y, .y = ab.z * ac.x - ab.x * ac.z, .z = ab.x * ac.y - ab.y * ac.x };
        const dot = n.x * a.x + n.y * a.y + n.z * a.z;
        std.debug.print("WIND: face={d} dot={d:.6}\n", .{ ff, dot });
        if (!std.math.signbit(sign_ref) and std.math.signbit(dot)) all_same = false;
        if (std.math.signbit(sign_ref) and !std.math.signbit(dot)) all_same = false;
    }
    try expect(all_same);
}

test "grid: pentas alias coverage log" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();

    const Ni: isize = @intCast(g.size);
    var face: usize = 0;
    var missing: usize = 0;
    var found: usize = 0;
    while (face < 20) : (face += 1) {
        var q: isize = -Ni;
        while (q <= Ni) : (q += 1) {
            var r: isize = -Ni;
            while (r <= Ni) : (r += 1) {
                if (!g.isValid(q, r) or !g.isPenta(q, r)) continue;
                const h = planet.Grid.hashCoord(q, r, face);
                const idx_opt = m.get(h);
                if (idx_opt) |ix| {
                    std.debug.print("FOUND PENTA ALIAS: (f={d}, q={d}, r={d}) -> idx={d}\n", .{ face, q, r, ix });
                    found += 1;
                } else {
                    std.debug.print("MISSING PENTA ALIAS: (f={d}, q={d}, r={d})\n", .{ face, q, r });
                    missing += 1;
                }
            }
        }
    }
    std.debug.print("SUMMARY: found={d} missing={d}\n", .{ found, missing });
}

test "grid: dump icosa faces and neighbors" {
    // Dump canonicalized icosa face vertex indices
    std.debug.print("ICOSAFACES (1-based vertex ids):\n", .{});
    var f: usize = 0;
    while (f < 20) : (f += 1) {
        const tri = gridmod.icosa_face_vertices[f];
        std.debug.print("F{d}: ({d},{d},{d})\n", .{ f, tri[0], tri[1], tri[2] });
    }

    // Build and dump neighbor faces
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    std.debug.print("NEIGHBORS (face -> nf0,nf1,nf2):\n", .{});
    f = 0;
    while (f < 20) : (f += 1) {
        const nf0 = g.faceNeighborFace(f, 0);
        const nf1 = g.faceNeighborFace(f, 1);
        const nf2 = g.faceNeighborFace(f, 2);
        std.debug.print("F{d}: {d},{d},{d}\n", .{ f, nf0, nf1, nf2 });
    }
}

test "grid: face_axes parity vs pentagon alias coverage" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();

    var m = try g.populateIndices();
    defer m.deinit();

    const Parity = enum { even, odd };

    const N: isize = @intCast(g.size);
    const PENTA6: [6][2]isize = .{ .{ -N, 0 }, .{ -N, N }, .{ 0, N }, .{ N, 0 }, .{ N, -N }, .{ 0, -N } };

    std.debug.print("FACE_AXES PARITY AND PENTA COVERAGE:\n", .{});
    var f: usize = 0;
    var complete_faces: [20]bool = .{ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false };
    while (f < 20) : (f += 1) {
        const perm = g.face_axes[f];
        var inv: u8 = 0;
        var i: u8 = 0;
        while (i < 3) : (i += 1) {
            var j: u8 = i + 1;
            while (j < 3) : (j += 1) {
                if (perm[i] > perm[j]) inv += 1;
            }
        }
        const parity: Parity = if ((inv & 1) == 0) .even else .odd;
        // Count missing pentagon aliases for this face
        var missing: usize = 0;
        var k: usize = 0;
        while (k < 6) : (k += 1) {
            const q = PENTA6[k][0];
            const r = PENTA6[k][1];
            // Only consider if recognized as a pentagon in canonical test (matches previous coverage test)
            if (!g.isValid(q, r) or !g.isPenta(q, r)) continue;
            const h = planet.Grid.hashCoord(q, r, f);
            if (m.get(h) == null) missing += 1;
        }
        complete_faces[f] = (missing == 0);
        std.debug.print("F{d}: perm=({d},{d},{d}) parity={s} missing={d}\n", .{ f, perm[0], perm[1], perm[2], if (parity == .odd) "odd" else "even", missing });
    }

    std.debug.print("COMPLETE_FACES: ", .{});
    f = 0;
    while (f < 20) : (f += 1) if (complete_faces[f]) std.debug.print("{d} ", .{f});
    std.debug.print("\n", .{});
}

test "grid: parity instrumentation at (f0,f5,q-3,r4,dir5)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    g.testSeamDebug(-3, 4, 0, 5);
    g.testSeamDebug(-3, 4, 5, 5);
}

test "grid: geom neighbor vs step for missing alias (f0,q-4,r4)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();

    var m = try g.populateIndices();
    defer m.deinit();

    // Missing alias we want to probe
    const f_m: usize = 0;
    const q_m: isize = -4;
    const r_m: isize = 4;

    // Compute 3D position of the missing alias
    const pos_m = g.faceCenterPosition(q_m, r_m, f_m);

    // Reference neighbor distance: use a known present pentagon tile and its nearest stepped neighbor
    const pref_q: isize = 0;
    const pref_r: isize = 4;
    const pref_f: usize = 5; // known present
    const pos_pref = g.faceCenterPosition(pref_q, pref_r, pref_f);
    var d_ref: f64 = 1e9;
    var dir: u8 = 0;
    while (dir < 6) : (dir += 1) {
        const nb = g.testStep(pref_q, pref_r, pref_f, dir);
        const pos_nb = g.faceCenterPosition(nb.q, nb.r, nb.face);
        const dx = pos_nb.x - pos_pref.x;
        const dy = pos_nb.y - pos_pref.y;
        const dz = pos_nb.z - pos_pref.z;
        const d = std.math.sqrt(dx * dx + dy * dy + dz * dz);
        if (d > 0 and d < d_ref) d_ref = d;
    }
    // Use a slightly generous threshold
    const thresh = d_ref * 1.15;

    // Iterate discovered tiles and find geometric neighbors of the missing alias
    var close_count: usize = 0;
    var i: usize = 0;
    while (i < g.tile_count) : (i += 1) {
        const c = g.coords[i];
        const pos_c = g.faceCenterPosition(c.q, c.r, c.face);
        const dx = pos_c.x - pos_m.x;
        const dy = pos_c.y - pos_m.y;
        const dz = pos_c.z - pos_m.z;
        const d = std.math.sqrt(dx * dx + dy * dy + dz * dz);
        if (d <= thresh) {
            close_count += 1;
            // Check if stepping from c reaches the missing alias
            var reaches = false;
            var d2: u8 = 0;
            while (d2 < 6) : (d2 += 1) {
                const nb = g.testStep(c.q, c.r, c.face, d2);
                if (nb.face == f_m and nb.q == q_m and nb.r == r_m) {
                    reaches = true;
                    break;
                }
            }
            std.debug.print("GEOM_NEI: c=(f{d},q{d},r{d}) d={d:.6} reaches_missing={any}\n", .{ c.face, c.q, c.r, d, reaches });
        }
    }
    std.debug.print("GEOM_NEI SUMMARY: found_close={d}, d_ref={d:.6}, thresh={d:.6}\n", .{ close_count, d_ref, thresh });
}

test "grid: two-hop around pentagon (f0,q-3,r4) to missing (f0,q-4,r4)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();

    // Missing alias and pentagon-adjacent start
    const f_m: usize = 0;
    const q_m: isize = -4;
    const r_m: isize = 4;
    const pos_m = g.faceCenterPosition(q_m, r_m, f_m);

    const f_p: usize = 0;
    const q_p: isize = -3;
    const r_p: isize = 4;
    std.debug.print("PENT_START: P=(f{d},q{d},r{d})  M=(f{d},q{d},r{d})\n", .{ f_p, q_p, r_p, f_m, q_m, r_m });

    // One-hop sweep from P
    var best_dir: u8 = 0;
    var best_d: f64 = 1e9;
    var dir: u8 = 0;
    while (dir < 6) : (dir += 1) {
        const nb1 = g.testStep(q_p, r_p, f_p, dir);
        const pos_nb1 = g.faceCenterPosition(nb1.q, nb1.r, nb1.face);
        const dx = pos_nb1.x - pos_m.x;
        const dy = pos_nb1.y - pos_m.y;
        const dz = pos_nb1.z - pos_m.z;
        const d = std.math.sqrt(dx * dx + dy * dy + dz * dz);
        const eq_missing = (nb1.face == f_m and nb1.q == q_m and nb1.r == r_m);
        std.debug.print("ONE_HOP: dir={d} -> nb1=(f{d},q{d},r{d}) d_to_M={d:.6} eq_missing={any}\n", .{ dir, nb1.face, nb1.q, nb1.r, d, eq_missing });
        if (d < best_d) {
            best_d = d;
            best_dir = dir;
        }
    }
    std.debug.print("ONE_HOP_BEST: dir={d} d={d:.6}\n", .{ best_dir, best_d });

    // Two-hop sweep from the best one-hop neighbor
    const nb_best = g.testStep(q_p, r_p, f_p, best_dir);
    var reached = false;
    var dir2: u8 = 0;
    while (dir2 < 6) : (dir2 += 1) {
        const nb2 = g.testStep(nb_best.q, nb_best.r, nb_best.face, dir2);
        const eq_missing2 = (nb2.face == f_m and nb2.q == q_m and nb2.r == r_m);
        const pos_nb2 = g.faceCenterPosition(nb2.q, nb2.r, nb2.face);
        const dx2 = pos_nb2.x - pos_m.x;
        const dy2 = pos_nb2.y - pos_m.y;
        const dz2 = pos_nb2.z - pos_m.z;
        const d2 = std.math.sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
        std.debug.print("TWO_HOP: dir1={d} dir2={d} -> nb2=(f{d},q{d},r{d}) d_to_M={d:.6} eq_missing={any}\n", .{ best_dir, dir2, nb2.face, nb2.q, nb2.r, d2, eq_missing2 });
        if (eq_missing2) reached = true;
    }
    std.debug.print("TWO_HOP_REACHED: {any}\n", .{reached});
}

test "grid: broad two-hop sweep to missing alias (f0,q-4,r4)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    var m = try g.populateIndices();
    defer m.deinit();

    // Missing alias target
    const f_m: usize = 0;
    const q_m: isize = -4;
    const r_m: isize = 4;
    const pos_m = g.faceCenterPosition(q_m, r_m, f_m);

    // Use a present pentagon to establish a neighbor-distance baseline
    const pref_q: isize = 0;
    const pref_r: isize = 4;
    const pref_f: usize = 5;
    const pos_pref = g.faceCenterPosition(pref_q, pref_r, pref_f);
    var d_ref: f64 = 1e9;
    var dir: u8 = 0;
    while (dir < 6) : (dir += 1) {
        const nb = g.testStep(pref_q, pref_r, pref_f, dir);
        const pos_nb = g.faceCenterPosition(nb.q, nb.r, nb.face);
        const dx = pos_nb.x - pos_pref.x;
        const dy = pos_nb.y - pos_pref.y;
        const dz = pos_nb.z - pos_pref.z;
        const d = std.math.sqrt(dx * dx + dy * dy + dz * dz);
        if (d > 0 and d < d_ref) d_ref = d;
    }
    const thresh = d_ref * 1.15;
    std.debug.print("SWEEP: d_ref={d:.6} thresh={d:.6}\n", .{ d_ref, thresh });

    // Gather geometric neighbors of the missing alias
    var candidates: usize = 0;
    var found_paths: usize = 0;
    var i: usize = 0;
    while (i < g.tile_count) : (i += 1) {
        const c = g.coords[i];
        const pos_c = g.faceCenterPosition(c.q, c.r, c.face);
        const dx = pos_c.x - pos_m.x;
        const dy = pos_c.y - pos_m.y;
        const dz = pos_c.z - pos_m.z;
        const d = std.math.sqrt(dx * dx + dy * dy + dz * dz);
        if (d > thresh) continue;
        candidates += 1;

        // Try all two-hop sequences from this candidate
        var d1: u8 = 0;
        while (d1 < 6) : (d1 += 1) {
            const mid = g.testStep(c.q, c.r, c.face, d1);
            var d2: u8 = 0;
            while (d2 < 6) : (d2 += 1) {
                const fin = g.testStep(mid.q, mid.r, mid.face, d2);
                const hit = (fin.face == f_m and fin.q == q_m and fin.r == r_m);
                if (hit) {
                    std.debug.print(
                        "PATH_FOUND: start=(f{d},q{d},r{d}) dir1={d} mid=(f{d},q{d},r{d}) dir2={d} -> (f{d},q{d},r{d})\n",
                        .{ c.face, c.q, c.r, d1, mid.face, mid.q, mid.r, d2, fin.face, fin.q, fin.r },
                    );
                    found_paths += 1;
                }
            }
        }
    }
    std.debug.print("SWEEP SUMMARY: candidates={d} found_paths={d}\n", .{ candidates, found_paths });
}

test "corner-turn row audit for face 0 @ (-4,0)" {
    var g = try planet.Grid.init(std.testing.allocator, 4);
    defer g.deinit();
    const f: usize = 0;
    const v = g.faceLocalVertexForCorner(f, -4, 0);
    const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
    var e: usize = 0;
    while (e < 3) : (e += 1) {
        const inc = (pairs[e][0] == v) or (pairs[e][1] == v);
        if (!inc) continue;
        g.debugCornerTurnEntry(f, v, e);
    }
}

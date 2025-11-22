const std = @import("std");
const Vec3 = struct { x: f64, y: f64, z: f64 };
const QR = struct { q: isize, r: isize };
const EdgeKey = struct { a: usize, b: usize };
fn edgeKey(a: usize, b: usize) EdgeKey {
    return .{ .a = @min(a, b), .b = @max(a, b) };
}
const UNSET = std.math.maxInt(usize);
pub const Coord = struct { q: isize, r: isize, face: usize };
pub const NeighborSet = struct { ids: [6]usize, count: u8 };
comptime {
    // 22 (face) + 21 (q) + 21 (r) = 64
    std.debug.assert((22 + 21 + 21) == 64);
}

// (debug helpers defined as Grid methods inside the struct)

inline fn legacyRemoved(comptime msg: []const u8) noreturn {
    @compileError(msg);
}

fn cross(a: Vec3, b: Vec3) Vec3 {
    return .{
        .x = a.y * b.z - a.z * b.y,
        .y = a.z * b.x - a.x * b.z,
        .z = a.x * b.y - a.y * b.x,
    };
}
fn dot(a: Vec3, b: Vec3) f64 {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

fn canonicalizeFacesComptime(comptime F_in: [20][3]usize) [20][3]usize {
    @setEvalBranchQuota(200000);
    var F = F_in; // mutable copy (1-based indices)

    // BFS over faces; take face 0 as canonical local order
    var seen = [_]bool{false} ** 20;
    var q: [20]usize = undefined;
    var head: usize = 0;
    var tail: usize = 0;
    seen[0] = true;
    q[tail] = 0;
    tail += 1;

    while (head < tail) {
        const f = q[head];
        head += 1;
        const tri = F[f];
        inline for (0..3) |ei| {
            const a = tri[ei];
            const b = tri[(ei + 1) % 3];
            // find neighbor that shares undirected edge {a,b}
            var nf: usize = 20;
            var ne: usize = 3;
            var dirSame: bool = false;
            var i: usize = 0;
            while (i < 20) : (i += 1) {
                if (i == f) continue;
                const gtest = F[i];
                inline for (0..3) |k| {
                    const x = gtest[k];
                    const y = gtest[(k + 1) % 3];
                    if (x == a and y == b) {
                        nf = i;
                        ne = k;
                        dirSame = true;
                        break;
                    }
                    if (x == b and y == a) {
                        nf = i;
                        ne = k;
                        dirSame = false;
                        break;
                    }
                }
                if (nf != 20) break;
            }
            if (nf == 20) continue;
            if (!seen[nf]) {
                if (dirSame) {
                    // neighbor lists a->b; we need b->a. Flip neighbor triangle.
                    const tmp = F[nf][1];
                    F[nf][1] = F[nf][2];
                    F[nf][2] = tmp;
                }
                // rotate neighbor to place consecutive (b,a) at (ei+1, ei) positions
                const g = F[nf];
                // find rotation such that (rot0, rot1) == (b,a)
                inline for (0..3) |rot| {
                    const x = g[(rot + 0) % 3];
                    const y = g[(rot + 1) % 3];
                    const z = g[(rot + 2) % 3];
                    if (x == b and y == a) {
                        var out = g;
                        out[(ei + 1) % 3] = x;
                        out[(ei + 2) % 3] = y;
                        out[(ei + 0) % 3] = z;
                        F[nf] = out;
                        break;
                    }
                }
                seen[nf] = true;
                q[tail] = nf;
                tail += 1;
            } else {
                // already seen: assert opposite direction on neighbor
                const g = F[nf];
                const x = g[ne];
                const y = g[(ne + 1) % 3];
                std.debug.assert(x == b and y == a);
            }
        }
    }

    // sanity: each undirected edge appears exactly twice
    inline for (0..20) |i| {
        const t = F[i];
        inline for (0..3) |k| {
            const u = t[k];
            const v = t[(k + 1) % 3];
            var count: usize = 0;
            inline for (0..20) |j| {
                const s = F[j];
                inline for (0..3) |m| {
                    const x = s[m];
                    const y = s[(m + 1) % 3];
                    if ((x == u and y == v) or (x == v and y == u)) count += 1;
                }
            }
            std.debug.assert(count == 2);
        }
    }
    return F;
}

pub const icosa_vertices = [_]Vec3{
    .{ .x = 0.850650808, .y = 0.525731112, .z = 0.000000000 }, // v0  = ( φ,  1, 0) / √(1+φ²)
    .{ .x = -0.850650808, .y = 0.525731112, .z = 0.000000000 }, // v1  = (-φ,  1, 0) / √(1+φ²)
    .{ .x = 0.850650808, .y = -0.525731112, .z = 0.000000000 }, // v2  = ( φ, -1, 0) / √(1+φ²)
    .{ .x = -0.850650808, .y = -0.525731112, .z = 0.000000000 }, // v3  = (-φ, -1, 0) / √(1+φ²)
    .{ .x = 0.525731112, .y = 0.000000000, .z = 0.850650808 }, // v4  = ( 1,  0, φ) / √(1+φ²)
    .{ .x = 0.525731112, .y = 0.000000000, .z = -0.850650808 }, // v5  = ( 1,  0,-φ) / √(1+φ²)
    .{ .x = -0.525731112, .y = 0.000000000, .z = 0.850650808 }, // v6  = (-1,  0, φ) / √(1+φ²)
    .{ .x = -0.525731112, .y = 0.000000000, .z = -0.850650808 }, // v7  = (-1,  0,-φ) / √(1+φ²)
    .{ .x = 0.000000000, .y = 0.850650808, .z = 0.525731112 }, // v8  = ( 0, φ,  1) / √(1+φ²)
    .{ .x = 0.000000000, .y = -0.850650808, .z = 0.525731112 }, // v9  = ( 0,-φ,  1) / √(1+φ²)
    .{ .x = 0.000000000, .y = 0.850650808, .z = -0.525731112 }, // v10 = ( 0, φ, -1) / √(1+φ²)
    .{ .x = 0.000000000, .y = -0.850650808, .z = -0.525731112 }, // v11 = ( 0,-φ, -1) / √(1+φ²)
};

// The 20 faces connecting the 12 vertices above (using 1-based indexing).
const icosa_face_vertices_raw = [_][3]usize{
    .{ 1, 9, 5 },   .{ 1, 6, 11 },  .{ 3, 5, 10 }, .{ 3, 12, 6 }, .{ 2, 7, 9 },
    .{ 2, 11, 8 },  .{ 4, 10, 7 },  .{ 4, 8, 12 }, .{ 1, 11, 9 }, .{ 2, 9, 11 },
    .{ 3, 10, 12 }, .{ 4, 10, 12 }, .{ 5, 3, 1 },  .{ 6, 1, 3 },  .{ 7, 2, 4 },
    .{ 8, 4, 2 },   .{ 9, 7, 5 },   .{ 10, 5, 7 }, .{ 11, 6, 8 }, .{ 12, 8, 6 },
};

// Computed at compile time: canonicalized local indices per face so shared edges align
pub const icosa_face_vertices = canonicalizeFacesComptime(icosa_face_vertices_raw);

// Face neighbors will be constructed at init from faces

// Canonical axial direction list (CW): (+q), (+r), (-q,+s), (-q), (-r), (+q,-s)
const DIRS = [_][2]isize{ .{ 1, 0 }, .{ 1, -1 }, .{ 0, -1 }, .{ -1, 0 }, .{ -1, 1 }, .{ 0, 1 } };

// Barycentric step deltas in LABEL space (A,B,C) for the 6 hex directions:
// d0: +A -B, d1: +A -C, d2: +B -C, d3: -A +B, d4: -A +C, d5: -B +C
const DIR_BARY_LABEL = [_][3]isize{
    .{ 1, -1, 0 },
    .{ 1, 0, -1 },
    .{ 0, 1, -1 },
    .{ -1, 1, 0 },
    .{ -1, 0, 1 },
    .{ 0, -1, 1 },
};

var debug_edge_logs: usize = 0;
var debug_nbr_logs: usize = 0;
var test_probe_done: bool = false;
var seam_verify_logged: bool = false;
// Lightweight instrumentation counters (debugging seam hops); no prints in hot paths
var g_seam_hop_calls: usize = 0;
var g_seam_project_calls: usize = 0;
var g_seam_iters: usize = 0;

pub const Grid = struct {
    allocator: std.mem.Allocator,
    size: usize,
    tile_count: usize,

    // Final topology tables
    coords: []Coord,
    neighbors: []NeighborSet,
    face_neighbors: [20][3][2]usize,
    // permutation mapping of barycentric coords across each face edge
    face_perm: [20][3][3]u8,

    // per-face axes mapping from face-local vertex order (0,1,2) to standard (u,v,w)
    face_axes: [20][3]u8,
    // Normalized (0-based) face label triples (A,B,C) after global winding+rotation normalization
    normalized_faces: [20][3]u8,
    // Guarded seam-perm writes
    _perm_writes: [20][3]u8 = std.mem.zeroes([20][3]u8),
    _seams_frozen: bool = false,
    // Guarded neighbor writes
    _nei_writes: [20][3]u8 = std.mem.zeroes([20][3]u8),
    _neighbors_frozen: bool = false,

    // Debug snapshots to guard against accidental seam rewrites after init (Debug mode only)
    _perm_snapshot: [20][3][3]u8 = std.mem.zeroes([20][3][3]u8),
    _ne_snapshot: [20][3][2]usize = std.mem.zeroes([20][3][2]usize),

    // Internal data structures for tracking shared indices during generation
    penta_indices: [12]usize,
    edge_indices: [20][3][]usize,
    // Deterministic ownership for boundary canonicalization
    edge_owner: [20][3]u8 = std.mem.zeroes([20][3]u8),
    vertex_owner: [12]u8 = std.mem.zeroes([12]u8),

    // Generated mesh data (flat arrays suitable for GPU upload)
    // vertices: xyz triplets in f32
    vertices: ?[]f32,
    // elevations per-vertex (f32)
    elevations: ?[]f32,
    // triangle indices (u32)
    indices: ?[]u32,

    // Public helper struct for edge mappings (used in tests)
    pub const EdgeMap = struct { nf: usize, ne: usize, rev: bool };

    // Sum-preserving rounding of barycentric floats scaled to 3N.
    // Floors each component, then distributes residual by largest fractional parts.

    // Geometric oracle: planar barycentric coordinates of 3D point P on face 'face_idx'
    pub fn baryFromPointOnFace(Px: f64, Py: f64, Pz: f64, face_idx: usize) [3]f64 {
        const tri = icosa_face_vertices[face_idx];
        const A = icosa_vertices[tri[0] - 1];
        const B = icosa_vertices[tri[1] - 1];
        const C = icosa_vertices[tri[2] - 1];
        const v0x = B.x - A.x;
        const v0y = B.y - A.y;
        const v0z = B.z - A.z;
        const v1x = C.x - A.x;
        const v1y = C.y - A.y;
        const v1z = C.z - A.z;
        const v2x = Px - A.x;
        const v2y = Py - A.y;
        const v2z = Pz - A.z;
        const d00 = v0x * v0x + v0y * v0y + v0z * v0z;
        const d01 = v0x * v1x + v0y * v1y + v0z * v1z;
        const d11 = v1x * v1x + v1y * v1y + v1z * v1z;
        const d20 = v2x * v0x + v2y * v0y + v2z * v0z;
        const d21 = v2x * v1x + v2y * v1y + v2z * v1z;
        const denom = d00 * d11 - d01 * d01;
        if (denom == 0) return .{ 1.0, 0.0, 0.0 };
        const v = (d11 * d20 - d01 * d21) / denom;
        const w = (d00 * d21 - d01 * d20) / denom;
        const u = 1.0 - v - w;
        return .{ u, v, w };
    }
    pub fn init(allocator: std.mem.Allocator, size: usize) !Grid {
        std.debug.assert(size >= 1);
        std.debug.assert(size <= (@as(usize, 1) << 20));
        var grid = Grid{
            .allocator = allocator,
            .size = size,
            .tile_count = 0,
            .coords = &[_]Coord{},
            .neighbors = &[_]NeighborSet{},
            .face_neighbors = std.mem.zeroes([20][3][2]usize),
            .face_perm = std.mem.zeroes([20][3][3]u8),
            .face_axes = std.mem.zeroes([20][3]u8),
            .normalized_faces = std.mem.zeroes([20][3]u8),
            ._nei_writes = std.mem.zeroes([20][3]u8),

            ._perm_snapshot = std.mem.zeroes([20][3][3]u8),
            ._ne_snapshot = std.mem.zeroes([20][3][2]usize),
            .penta_indices = undefined,
            .edge_indices = std.mem.zeroes([20][3][]usize),
            .vertices = null,
            .elevations = null,
            .indices = null,
        };
        // Initialize face_perm with sentinel 255 and clear write guards
        for (0..20) |f| {
            for (0..3) |e| {
                grid.face_perm[f][e] = .{ 255, 255, 255 };
                grid._perm_writes[f][e] = 0;
                grid._nei_writes[f][e] = 0;
            }
        }
        grid._seams_frozen = false;
        grid._neighbors_frozen = false;

        // Initialize penta sentinel and allocate edge indices arrays
        grid.penta_indices = .{UNSET} ** 12;
        for (0..20) |face| {
            for (0..3) |edge| {
                const edge_len: usize = if (size > 0) size - 1 else 0;
                grid.edge_indices[face][edge] = try allocator.alloc(usize, edge_len);
                if (edge_len > 0) @memset(grid.edge_indices[face][edge], UNSET);
            }
        }

        // Build label data → neighbors, then freeze neighbors
        grid.rebuildFaceTablesFromIcosa(); // neighbors built and frozen inside
        // Propagate axes across all faces from labels (seed face 0 = identity)
        grid.propagateAxesFromLabels();
        // Assert axes are fully set
        for (0..20) |fi| {
            const ax = grid.face_axes[fi];
            std.debug.assert(ax[0] < 3 and ax[1] < 3 and ax[2] < 3);
        }
        // Build perms from axes, then freeze seams
        grid.recomputeSeamPermsFromAxes();
        grid.freezeSeams();
        // Sanity: seam tables should be consistent

        // Lightweight: single-pass propagation, no repairs or verification
        grid.propagateAxesNoAssert();
        grid.recomputeSeamPermsFromAxes(); // now that face_axes exist, compute perms in actual basis
        return grid;
    }

    fn freezeSeamsForDebug(self: *Grid) void {
        if (@import("builtin").mode != .Debug) return;
        self._perm_snapshot = self.face_perm;
        self._ne_snapshot = self.face_neighbors;
    }

    fn dedupSorted(ids: []usize) usize {
        if (ids.len == 0) return 0;
        var w: usize = 1;
        var i: usize = 1;
        while (i < ids.len) : (i += 1) {
            if (ids[i] != ids[i - 1]) {
                ids[w] = ids[i];
                w += 1;
            }
        }
        if (w >= 2 and ids[0] == ids[w - 1]) w -= 1;
        return w;
    }

    pub fn deinit(self: *Grid) void {
        if (self.coords.len > 0) self.allocator.free(self.coords);
        if (self.neighbors.len > 0) self.allocator.free(self.neighbors);
        for (0..20) |face| {
            for (0..3) |edge| {
                self.allocator.free(self.edge_indices[face][edge]);
            }
        }

        if (self.vertices) |v| self.allocator.free(v);
        if (self.elevations) |e| self.allocator.free(e);
        if (self.indices) |i| self.allocator.free(i);
    }

    fn buildPermForEdge(self: *const Grid, f: usize, e: usize) [3]u8 {
        // Build CURRENT face-local -> NEIGHBOR face-local permutation for seam (f,e) from geometry.
        const pairs = [_][2]usize{ .{ 1, 2 }, .{ 2, 0 }, .{ 0, 1 } };
        const nf = self.face_neighbors[f][e][0] - 1;
        const cur = icosa_face_vertices[f];
        const nfv = icosa_face_vertices[nf];
        const a_local = pairs[e][0];
        const b_local = pairs[e][1];
        const ga = cur[a_local];
        const gb = cur[b_local];
        var ia: usize = 3;
        var ib: usize = 3;
        var i: usize = 0;
        while (i < 3) : (i += 1) {
            if (nfv[i] == ga) ia = i;
            if (nfv[i] == gb) ib = i;
        }
        if (ia == 3 or ib == 3) {
            if (@import("builtin").is_test) {
                std.debug.print("buildPermForEdge FAIL: f={d} e={d} nf={d} ga={d} gb={d} cur=({d},{d},{d}) nfv=({d},{d},{d})\n", .{
                    f, e, nf, ga, gb, cur[0], cur[1], cur[2], nfv[0], nfv[1], nfv[2],
                });
            }
            @panic("buildPermForEdge: shared vertices not found");
        }
        const c_local: usize = 3 - a_local - b_local;
        const ic: usize = 3 - ia - ib;
        // Build dest->src mapping for applyPerm3: uvw_nf[j] = uvw_src[perm[j]]
        var perm_fwd: [3]u8 = .{ 0, 0, 0 };
        perm_fwd[@intCast(ia)] = @intCast(a_local);
        perm_fwd[@intCast(ib)] = @intCast(b_local);
        perm_fwd[@intCast(ic)] = @intCast(c_local);
        return perm_fwd;
    }

    fn rebuildFaceTablesFromIcosa(self: *Grid) void {
        // Normalize faces directly from canonicalized icosa faces (0-based labels)
        var faces0: [20][3]u8 = undefined;
        for (0..20) |i| {
            faces0[i][0] = @intCast(icosa_face_vertices[i][0] - 1);
            faces0[i][1] = @intCast(icosa_face_vertices[i][1] - 1);
            faces0[i][2] = @intCast(icosa_face_vertices[i][2] - 1);
        }
        self.normalized_faces = faces0;
        // Normalize label rotations/winding via BFS to enforce oriented-reversed seams
        self.normalizeFacesByBfs();
        // Rebuild directed reversed-edge neighbors from labels
        for (0..20) |i| {
            for (0..3) |e| {
                self.face_neighbors[i][e] = .{ 0, 0 };
            }
        }
        self.rebuildNeighborsFromNormalizedLabels();
        // Freeze neighbors before any seam building
        self.freezeNeighbors();
        // Build deterministic edge owners now that neighbors are frozen
        self.buildEdgeOwners();
        self.buildVertexOwners();
        if (@import("builtin").is_test and TEST_VERBOSE) {
            const f0 = self.normalized_faces[0];
            const f16_lbl = self.normalized_faces[16];
            std.debug.print("NORM_FACES: f0=({d},{d},{d}) f16=({d},{d},{d})\n", .{ f0[0], f0[1], f0[2], f16_lbl[0], f16_lbl[1], f16_lbl[2] });
            const nf_e0 = self.face_neighbors[0][0][0] - 1;
            const ne_e0 = self.face_neighbors[0][0][1] - 1;
            std.debug.print("NEI_F0_E0: nf={d} ne={d}\n", .{ nf_e0, ne_e0 });
        }
    }

    // Rebuild neighbor mapping from the current normalized_faces using only labels (directed reversed-edge).
    fn rebuildNeighborsFromNormalizedLabels(self: *Grid) void {
        for (0..20) |i| {
            for (0..3) |e| self.face_neighbors[i][e] = .{ 0, 0 };
        }
        var fi: usize = 0;
        while (fi < 20) : (fi += 1) {
            const Pf = self.normalized_faces[fi];
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                if (self._nei_writes[fi][e] != 0) continue; // already set (by reciprocal earlier)
                const a: u8 = Pf[(e + 1) % 3];
                const b: u8 = Pf[(e + 2) % 3];
                const nf = self.findFaceByTwoLabels(a, b, fi);
                const Pn = self.normalized_faces[nf];
                const ne = findReversedOnNF(Pn, a, b) orelse @panic("oriented nf edge not found");
                self.setNeighbor(fi, e, nf, ne, "oriented-builder");
                self.setNeighbor(nf, ne, fi, e, "oriented-builder-recip");
            }
        }
        // Reciprocity
        for (0..20) |ff| {
            for (0..3) |ee| {
                const nfa = self.face_neighbors[ff][ee][0];
                const nea = self.face_neighbors[ff][ee][1];
                if (nfa == 0 or nea == 0) continue; // skip unset
                const nfo = nfa - 1;
                const neo = nea - 1;
                const back = self.face_neighbors[nfo][neo];
                std.debug.assert(back[0] == ff + 1 and back[1] == ee + 1);
            }
        }
    }

    fn recomputeSeamPermsFromAxes(self: *Grid) void {
        // Build-only: neighbors are frozen; axes are propagated.
        // Write both directions per seam once; skip if reciprocal already set.
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                if (self.face_perm[f][e][0] != 255) continue;
                const nf0 = self.face_neighbors[f][e][0];
                const ne0 = self.face_neighbors[f][e][1];
                if (nf0 == 0 or ne0 == 0) continue; // skip unset slots defensively
                const nf = nf0 - 1;
                const ne = ne0 - 1;
                // Re-derive neighbor edge from labels (reversed oriented edge) and assert
                const Pf_lbl = self.normalized_faces[f];
                const Pnf_lbl = self.normalized_faces[nf];
                const pairsF = [3][2]u8{ .{ 1, 2 }, .{ 2, 0 }, .{ 0, 1 } };
                const ix_f: u8 = pairsF[e][0];
                const iy_f: u8 = pairsF[e][1];
                const ne_expected = edgeIndexForPairOrientedLocalLbl(Pnf_lbl, Pf_lbl[iy_f], Pf_lbl[ix_f]);
                std.debug.assert(ne_expected < 3 and ne == ne_expected);
                // Compute p_actual directly from axes (dst actual on nf -> src actual on f)
                const p_actual = expectedPermFromAxes(self.face_axes[f], self.face_axes[nf]);
                self.setSeamPerm(f, e, p_actual, "builder f,e");
                self.setSeamPerm(nf, ne, invertPerm3(p_actual), "builder nf,ne");
            }
        }
        self.freezeSeams();
    }

    inline fn applyPerm3(x: [3]isize, p: [3]u8) [3]isize {
        return .{ x[p[0]], x[p[1]], x[p[2]] };
    }

    inline fn indexOf3u8(a: [3]u8, needle: u8) u8 {
        if (a[0] == needle) return 0;
        if (a[1] == needle) return 1;
        return 2;
    }
    inline fn invertPerm3(p: [3]u8) [3]u8 {
        // Harden against invalid inputs during early build phases
        if (p[0] > 2 or p[1] > 2 or p[2] > 2) return .{ 0, 1, 2 };
        var inv: [3]u8 = .{ 0, 0, 0 };
        inv[p[0]] = 0;
        inv[p[1]] = 1;
        inv[p[2]] = 2;
        return inv;
    }
    inline fn isValidPerm3(p: [3]u8) bool {
        if (p[0] > 2 or p[1] > 2 or p[2] > 2) return false;
        return (p[0] != p[1] and p[0] != p[2] and p[1] != p[2]);
    }
    inline fn setNeighbor(self: *Grid, f: usize, e: usize, nf: usize, ne: usize, reason: []const u8) void {
        std.debug.assert(!self._neighbors_frozen);
        std.debug.assert(f < 20 and e < 3 and nf < 20 and ne < 3);
        std.debug.assert(self._nei_writes[f][e] == 0);
        self.face_neighbors[f][e][0] = nf + 1;
        self.face_neighbors[f][e][1] = ne + 1;
        self._nei_writes[f][e] = 1;
        if (@import("builtin").is_test and TEST_VERBOSE) {
            std.debug.print("SETNEI f={d} e={d} -> nf={d} ne={d} from {s}\n", .{ f, e, nf, ne, reason });
        }
    }
    fn freezeNeighbors(self: *Grid) void {
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                std.debug.assert(self._nei_writes[f][e] == 1);
                const nei = self.face_neighbors[f][e];
                std.debug.assert(nei[0] != 0 and nei[1] != 0);
                const nf = nei[0] - 1;
                const ne = nei[1] - 1;
                std.debug.assert(nf < 20 and ne < 3);
            }
        }
        self._neighbors_frozen = true;
    }
    inline fn setSeamPerm(self: *Grid, f: usize, e: usize, p: [3]u8, reason: []const u8) void {
        std.debug.assert(!self._seams_frozen);
        std.debug.assert(f < 20 and e < 3);
        std.debug.assert(self._perm_writes[f][e] == 0);
        std.debug.assert(isValidPerm3(p));
        self.face_perm[f][e] = p;
        self._perm_writes[f][e] = 1;
        if (@import("builtin").is_test and TEST_VERBOSE) {
            std.debug.print("SETPERM f={d} e={d} p=({d},{d},{d}) from {s}\n", .{ f, e, p[0], p[1], p[2], reason });
        }
    }
    fn freezeSeams(self: *Grid) void {
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nei = self.face_neighbors[f][e];
                if (nei[0] == 0 or nei[1] == 0) continue; // unused seam slot; ignore
                std.debug.assert(self._perm_writes[f][e] == 1);
                const p = self.face_perm[f][e];
                std.debug.assert(isValidPerm3(p));
            }
        }
        self._seams_frozen = true;
    }
    fn buildEdgeOwners(self: *Grid) void {
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nei = self.face_neighbors[f][e];
                std.debug.assert(nei[0] != 0 and nei[1] != 0);
                const nf = nei[0] - 1;
                const ne = nei[1] - 1;
                const owner: u8 = @intCast(if (f < nf) f else nf);
                self.edge_owner[f][e] = owner;
                self.edge_owner[nf][ne] = owner;
            }
        }
    }
    fn buildVertexOwners(self: *Grid) void {
        // init to 255 (unset)
        var i: usize = 0;
        while (i < 12) : (i += 1) self.vertex_owner[i] = 255;
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            const P = self.normalized_faces[f];
            inline for (0..3) |k| {
                const v = P[k];
                const cur = self.vertex_owner[v];
                if (cur == 255 or f < cur) {
                    self.vertex_owner[v] = @intCast(f);
                }
            }
        }
        i = 0;
        while (i < 12) : (i += 1) {
            std.debug.assert(self.vertex_owner[i] != 255);
        }
    }
    inline fn faceHasLabel(self: *const Grid, f: usize, v: u8) bool {
        const P = self.normalized_faces[f];
        return (P[0] == v or P[1] == v or P[2] == v);
    }

    inline fn rotateTripleLbl(t: [3]u8, r: u2) [3]u8 {
        return switch (r) {
            0 => t,
            1 => .{ t[1], t[2], t[0] },
            2 => .{ t[2], t[0], t[1] },
            else => t,
        };
    }
    inline fn buildPermByLabelEquality(self: *Grid, f: usize, nf: usize) [3]u8 {
        // Build dst->src perm in ACTUAL basis by matching the two shared labels first,
        // then fill the remaining (opposite) axis deterministically.
        const Pf = self.normalized_faces[f];
        const Pnf = self.normalized_faces[nf];
        const Af = self.face_axes[f]; // label -> actual
        const invAf = invertPerm3(Af); // actual -> label
        const invAnf = invertPerm3(self.face_axes[nf]); // actual -> label
        var p_actual: [3]u8 = .{ 255, 255, 255 }; // dst (nf actual j) -> src (f actual i)
        var used_i: [3]bool = .{ false, false, false };
        // Map shared labels
        var li: u8 = 0;
        while (li < 3) : (li += 1) {
            const L = Pf[li];
            // Is this label present on nf?
            var present_on_nf = false;
            var j: u8 = 0;
            while (j < 3) : (j += 1) {
                if (Pnf[j] == L) {
                    present_on_nf = true;
                    break;
                }
            }
            if (!present_on_nf) continue; // skip the opposite (non-shared) label
            // Find actual axis j on nf carrying L
            var j_actual: u8 = 0;
            while (j_actual < 3) : (j_actual += 1) {
                if (Pnf[invAnf[j_actual]] == L) break;
            }
            // Find actual axis i on f carrying L
            var i_actual: u8 = 0;
            while (i_actual < 3) : (i_actual += 1) {
                if (Pf[invAf[i_actual]] == L) break;
            }
            p_actual[j_actual] = i_actual;
            used_i[i_actual] = true;
        }
        // Fill the remaining (opposite) axis
        var j_rem: u8 = 0;
        while (j_rem < 3 and p_actual[j_rem] != 255) : (j_rem += 1) {}
        var i_rem: u8 = 0;
        while (i_rem < 3 and used_i[i_rem]) : (i_rem += 1) {}
        if (j_rem < 3 and i_rem < 3) {
            p_actual[j_rem] = i_rem;
        }
        // Replace any leftover sentinels with a safe default to avoid OOB; should not happen
        var jj: u8 = 0;
        while (jj < 3) : (jj += 1) {
            if (p_actual[jj] == 255) p_actual[jj] = 0;
        }
        return p_actual;
    }

    inline fn edgeIndexForPairOrientedLocalLbl(t: [3]u8, a: u8, b: u8) u8 {
        if (t[1] == a and t[2] == b) return 0;
        if (t[2] == a and t[0] == b) return 1;
        if (t[0] == a and t[1] == b) return 2;
        return 3;
    }

    inline fn indexOfLabelInFaceU8(P: [3]u8, label: u8) u8 {
        if (P[0] == label) return 0;
        if (P[1] == label) return 1;
        if (P[2] == label) return 2;
        return 3;
    }
    inline fn indexOfLabelInFaceLbl(P: [3]usize, label: usize) u8 {
        var i: u8 = 0;
        while (i < 3) : (i += 1) {
            if (P[i] == label) return i;
        }
        return 3;
    }
    inline fn rot3(lbl: [3]u8, r: u8) [3]u8 {
        return .{ lbl[(0 + r) % 3], lbl[(1 + r) % 3], lbl[(2 + r) % 3] };
    }
    inline fn findReversedOnNF(Pn: [3]u8, a: u8, b: u8) ?u8 {
        var k: u8 = 0;
        while (k < 3) : (k += 1) {
            if (Pn[(k + 1) % 3] == b and Pn[(k + 2) % 3] == a) return k;
        }
        return null;
    }
    fn findFaceByTwoLabels(self: *const Grid, a: u8, b: u8, exclude: usize) usize {
        var k: usize = 0;
        while (k < 20) : (k += 1) {
            if (k == exclude) continue;
            const t = self.normalized_faces[k];
            var seen: u8 = 0;
            if (t[0] == a or t[0] == b) seen += 1;
            if (t[1] == a or t[1] == b) seen += 1;
            if (t[2] == a or t[2] == b) seen += 1;
            if (seen == 2) return k;
        }
        return exclude; // sentinel; caller should assert not equal
    }
    fn normalizeFacesByBfs(self: *Grid) void {
        // Enforce: for every seam (f,e)->nf, the oriented reversed edge exists on nf labels.
        var seen: [20]bool = .{false} ** 20;
        var queue: [20]u8 = undefined;
        var head: usize = 0;
        var tail: usize = 0;
        seen[0] = true;
        queue[tail] = 0;
        tail += 1;
        while (head < tail) : (head += 1) {
            const f = queue[head];
            const Pf = self.normalized_faces[f];
            var e: u8 = 0;
            while (e < 3) : (e += 1) {
                const a: u8 = Pf[(e + 1) % 3];
                const b: u8 = Pf[(e + 2) % 3];
                const nf = self.findFaceByTwoLabels(a, b, f);
                std.debug.assert(nf != f);
                if (!seen[nf]) {
                    // Rotate (or flip+rotate) nf so that (b,a) is an oriented edge
                    const Pn = self.normalized_faces[nf];
                    var ok = false;
                    var r: u8 = 0;
                    while (r < 3 and !ok) : (r += 1) {
                        const Pr = rot3(Pn, r);
                        if (findReversedOnNF(Pr, a, b) != null) {
                            self.normalized_faces[nf] = Pr;
                            ok = true;
                        }
                    }
                    if (!ok) {
                        const flip = .{ Pn[0], Pn[2], Pn[1] };
                        var rr: u8 = 0;
                        while (rr < 3 and !ok) : (rr += 1) {
                            const Pr = rot3(flip, rr);
                            if (findReversedOnNF(Pr, a, b) != null) {
                                self.normalized_faces[nf] = Pr;
                                ok = true;
                            }
                        }
                    }
                    std.debug.assert(ok);
                    seen[nf] = true;
                    queue[tail] = @intCast(nf);
                    tail += 1;
                } else {
                    const Pn = self.normalized_faces[nf];
                    std.debug.assert(findReversedOnNF(Pn, a, b) != null);
                }
            }
        }
    }
    inline fn buildPLblByEquality(Pf: [3]usize, Pnf: [3]usize) [3]u8 {
        // dst (nf label j) -> src (f label i) by label equality; fill remaining by elimination
        var p_lbl: [3]u8 = .{ 255, 255, 255 };
        var used_i: [3]bool = .{ false, false, false };
        var j: u8 = 0;
        while (j < 3) : (j += 1) {
            const L = Pnf[j];
            const i: u8 = indexOfLabelInFaceLbl(Pf, L);
            if (i < 3) {
                p_lbl[j] = i;
                used_i[i] = true;
            }
        }
        // Fill any remaining slots by the unused source index to complete a permutation
        var i_rem: u8 = 0;
        while (i_rem < 3 and used_i[i_rem]) : (i_rem += 1) {}
        var jj: u8 = 0;
        while (jj < 3) : (jj += 1) {
            if (p_lbl[jj] == 255) {
                p_lbl[jj] = if (i_rem < 3) i_rem else 0;
                used_i[p_lbl[jj]] = true;
            }
        }
        // validate permutation
        std.debug.assert(p_lbl[0] < 3 and p_lbl[1] < 3 and p_lbl[2] < 3);
        std.debug.assert(p_lbl[0] != p_lbl[1] and p_lbl[0] != p_lbl[2] and p_lbl[1] != p_lbl[2]);
        return p_lbl;
    }
    inline fn buildPLblFromReversedEdge(Pf: [3]usize, Pnf: [3]usize, e: u8) [3]u8 {
        // Oriented endpoints on f: e=0:(B,C), e=1:(C,A), e=2:(A,B)
        const pairsF = [3][2]u8{ .{ 1, 2 }, .{ 2, 0 }, .{ 0, 1 } };
        const ix_f: u8 = pairsF[e][0];
        const iy_f: u8 = pairsF[e][1];
        const x = Pf[ix_f];
        const y = Pf[iy_f];
        // On nf, match reversed oriented edge: (y,x)
        const iy_nf: u8 = indexOfLabelInFaceLbl(Pnf, y);
        const ix_nf: u8 = indexOfLabelInFaceLbl(Pnf, x);
        // Opposite component indices
        const ic_f: u8 = 3 - ix_f - iy_f;
        const ic_nf: u8 = 3 - ix_nf - iy_nf;
        var p_lbl: [3]u8 = .{ 0, 0, 0 }; // dst->src in STD basis
        // Swap edge endpoints (reversed)
        p_lbl[iy_nf] = ix_f; // nf(y) comes from f(x)
        p_lbl[ix_nf] = iy_f; // nf(x) comes from f(y)
        // Opposite ↔ opposite
        p_lbl[ic_nf] = ic_f;
        return p_lbl;
    }
    /// projectSeamIntoWindow
    ///
    /// Preconditions:
    /// - v is in label basis on the neighbor face
    /// - v[0] + v[1] + v[2] == R
    ///
    /// Seam aliases can produce tuples of the form (-a, -b, R + a + b),
    /// i.e. exactly two negatives, one coord > R, with sum fixed.
    /// We map that canonically to (0, 0, R) (up to permutation) by:
    ///   v_min0 += a; v_min1 += b; v_max -= (a + b);
    /// which preserves sum = R and lands in [0..R]^3.
    inline fn projectSeamIntoWindow(v: *[3]isize, R: isize) void {
        // Expect sum invariant
        if (@import("builtin").is_test) {
            std.debug.assert(v.*[0] + v.*[1] + v.*[2] == R);
        }
        // Early-exit if already within [0..R]^3
        if (v.*[0] >= 0 and v.*[0] <= R and v.*[1] >= 0 and v.*[1] <= R and v.*[2] >= 0 and v.*[2] <= R) {
            return;
        }
        // Sort indices by value ascending
        var idx0: usize = 0;
        var idx1: usize = 1;
        var idx2: usize = 2;
        if (v.*[idx0] > v.*[idx1]) {
            const tmp = idx0;
            idx0 = idx1;
            idx1 = tmp;
        }
        if (v.*[idx1] > v.*[idx2]) {
            const tmp = idx1;
            idx1 = idx2;
            idx2 = tmp;
        }
        if (v.*[idx0] > v.*[idx1]) {
            const tmp = idx0;
            idx0 = idx1;
            idx1 = tmp;
        }
        const imin0 = idx0; // most negative
        const imin1 = idx1; // second most negative
        const imax = idx2; // largest
        // Only handle canonical seam case (-,-,>R)
        if (!(v.*[imin0] < 0 and v.*[imin1] < 0 and v.*[imax] > R)) {
            return;
        }
        const a0 = -v.*[imin0];
        const a1 = -v.*[imin1];
        v.*[imin0] += a0; // -> 0
        v.*[imin1] += a1; // -> 0
        v.*[imax] -= (a0 + a1); // -> R
        if (@import("builtin").is_test) {
            std.debug.assert(v.*[0] >= 0 and v.*[0] <= R and v.*[1] >= 0 and v.*[1] <= R and v.*[2] >= 0 and v.*[2] <= R);
            std.debug.assert(v.*[0] + v.*[1] + v.*[2] == R);
        }
    }
    // Deterministic label-space seam hop with single-fold using endpoint donors only.
    inline fn seamHopLabelSpace(
        self: *const Grid,
        v_f_lab: [3]isize,
        dir: u8,
        f: usize,
        e: usize,
        nf: usize,
        ne: usize,
        p_lbl: [3]u8,
        R: isize,
    ) [3]isize {
        _ = .{ self, v_f_lab, dir, f, e, nf, ne, p_lbl, R };
        legacyRemoved("seamHopLabelSpace removed; seams are prebaked offline.");
        unreachable;
    }

    inline fn toBary(N: isize, q: isize, r: isize) [3]isize {
        const u = q + N;
        const v = r + N;
        const w = N - q - r;
        return .{ u, v, w };
    }

    // Public wrapper to allow external callers (e.g. prebaked builder) to use the same
    // barycentric conversion without exposing internal helpers broadly.
    pub inline fn toBary_public(N: isize, q: isize, r: isize) [3]isize {
        return toBary(N, q, r);
    }

    inline fn fromBary(N: isize, uvw: [3]isize) QR {
        return QR{ .q = uvw[0] - N, .r = uvw[1] - N };
    }

    // Very-near index resolver by 3D proximity on a small set of faces (additive-only usage).

    inline fn edgeIndexFromOppositeComponent(cn: usize) usize {
        // Low-bound rule: component index equals edge index:
        // cn=0 (A<0) -> e=0 (B,C)
        // cn=1 (B<0) -> e=1 (C,A)
        // cn=2 (C<0) -> e=2 (A,B)
        return cn;
    }

    inline fn isFaceHexStrict(_: *const Grid, uvw_face: [3]isize, N: isize) bool {
        if (uvw_face[0] < 0 or uvw_face[1] < 0 or uvw_face[2] < 0) return false;
        if (uvw_face[0] > 2 * N or uvw_face[1] > 2 * N or uvw_face[2] > 2 * N) return false;
        if (uvw_face[0] - uvw_face[2] > N) return false;
        if (uvw_face[1] - uvw_face[0] > N) return false;
        if (uvw_face[2] - uvw_face[1] > N) return false;
        return true;
    }

    inline fn oppositeComponentFromEdge(edge: usize) usize {
        // edges: 0:(0,1)->opp=2, 1:(1,2)->opp=0, 2:(2,0)->opp=1
        return switch (edge) {
            0 => 2,
            1 => 0,
            2 => 1,
            else => 0,
        };
    }

    pub fn canonicalizeGlobalAlias(self: *const Grid, face_in: usize, q_in: isize, r_in: isize) Coord {
        _ = .{ self, face_in, q_in, r_in };
        legacyRemoved("canonicalizeGlobalAlias removed; use prebaked TileId graph.");
        unreachable;
    }

    inline fn canonicalKey(self: *const Grid, face: usize, q: isize, r: isize) u64 {
        const c = self.canonicalizeGlobalAlias(face, q, r);
        return hashCoord(c.q, c.r, c.face);
    }

    pub inline fn canonicalCoord(self: *const Grid, face: usize, q: isize, r: isize) Coord {
        _ = .{ self, face, q, r };
        legacyRemoved("canonicalCoord removed; use prebaked TileId graph.");
        unreachable;
    }

    inline fn coordEq(_: *const Grid, a: Coord, b: Coord) bool {
        return a.face == b.face and a.q == b.q and a.r == b.r;
    }

    // Assert seam tables and stepper are mutually consistent without BFS.
    pub fn assertSeamTablesConsistent(self: *Grid) !void {
        const N: isize = @intCast(self.size);
        // Helper: pick a safe interior start (prefer 2, else 1, else 0)
        const base: isize = if (N >= 3) 2 else if (N >= 2) 1 else 0;
        const pairs = [_][2]usize{ .{ 1, 2 }, .{ 2, 0 }, .{ 0, 1 } };
        // Pass 1: stepper round-trip per (f,d) across first encountered seam
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            // Skip faces whose axes are not yet resolved
            if (!self.axesSet(f)) {
                continue;
            }
            var d: u8 = 0;
            while (d < 6) : (d += 1) {
                var start = Coord{ .face = f, .q = base, .r = base };
                var last = start;
                var k: usize = 0;
                var crossed = false;
                // Walk up to 2N+4 steps to ensure we reach a seam
                const max_steps: usize = @intCast(2 * N + 4);
                while (k < max_steps) : (k += 1) {
                    const nxt = self.stepPermuteExact(last.face, last.q, last.r, d);
                    if (nxt.face != last.face) {
                        start = last; // tile before crossing
                        crossed = true;
                        break;
                    }
                    last = nxt;
                }
                if (!crossed) continue; // nothing to check if we never hit a seam (tiny grids)
                // Compute violated edge on start face for dir d
                const uvw0 = self.toBaryFace(start.face, start.q, start.r);
                const duvw = self.dirCubeStep(start.face, d);
                const uvw1 = .{ uvw0[0] + duvw[0], uvw0[1] + duvw[1], uvw0[2] + duvw[2] };
                const cn_opt = self.pickViolatedEdgePost(start.face, uvw1, duvw);
                if (cn_opt == null) continue;
                const e = edgeIndexFromOppositeComponent(cn_opt.?);
                const nf = self.face_neighbors[start.face][e][0] - 1;
                // Step forward once across the seam
                const fwd = self.stepPermuteExact(start.face, start.q, start.r, d);
                std.debug.assert(fwd.face == nf);
                // Reverse dir: match -duvw transported to nf
                const perm = self.face_perm[start.face][e];
                const duvw_nf = applyPerm3(duvw, perm);
                const rdu = -duvw_nf[0];
                const rdv = -duvw_nf[1];
                const rdw = -duvw_nf[2];
                var rev_d: u8 = 255;
                var kk: u8 = 0;
                while (kk < 6) : (kk += 1) {
                    const cd = self.dirCubeStep(nf, kk);
                    if (cd[0] == rdu and cd[1] == rdv and cd[2] == rdw) {
                        rev_d = kk;
                        break;
                    }
                }
                std.debug.assert(rev_d != 255);
                const back = self.stepPermuteExact(fwd.face, fwd.q, fwd.r, rev_d);
                const back_can = self.canonicalCoord(back.face, back.q, back.r);
                const start_can = self.canonicalCoord(start.face, start.q, start.r);
                std.debug.assert(self.coordEq(back_can, start_can));
                // Neighbor-edge index and inverse perm consistency
                const perm_fwd = self.face_perm[start.face][e];
                // Compute neighbor edge index from endpoints via perm_fwd (find j where perm[j]==i)
                const a_i = pairs[e][0];
                const b_i = pairs[e][1];
                const ia = indexOf3u8(perm_fwd, @intCast(a_i));
                const ib = indexOf3u8(perm_fwd, @intCast(b_i));
                const ne_tab = self.face_neighbors[start.face][e][1] - 1;
                // Optional debug if derived-by-perm ne differs from stored ne
                const ne_calc = edgeIndexForPair(ia, ib);
                if (ne_calc != ne_tab and @import("builtin").is_test) {
                    std.debug.print("SEAM MISMATCH f={d} e={d} nf={d} ne_calc={d} ne_tab={d} ia={d} ib={d}\n", .{ start.face, e, nf, ne_calc, ne_tab, ia, ib });
                }
                // Check reverse neighbor mapping
                std.debug.assert(self.face_neighbors[nf][ne_tab][0] - 1 == start.face);
                std.debug.assert(self.face_neighbors[nf][ne_tab][1] - 1 == e);
                // Check inverse permutation
                const perm_bwd = self.face_perm[nf][ne_tab];
                std.debug.assert(perm_bwd[perm_fwd[0]] == 0 and perm_bwd[perm_fwd[1]] == 1 and perm_bwd[perm_fwd[2]] == 2);
            }
        }
        // Pass 2: geometric ground-truth pairing using ICOSAFACES order
        f = 0;
        while (f < 20) : (f += 1) {
            const tri = icosa_face_vertices[f];
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                const ne = self.face_neighbors[f][e][1] - 1;
                std.debug.assert(nf < 20 and ne < 3);
                const a_local = pairs[e][0];
                const b_local = pairs[e][1];
                const va = tri[a_local];
                const vb = tri[b_local];
                // find neighbor edge whose ordered endpoints are (vb,va)
                const ntri = icosa_face_vertices[nf];
                var ne_geom: usize = 3;
                var ee: usize = 0;
                while (ee < 3) : (ee += 1) {
                    const na = pairs[ee][0];
                    const nb = pairs[ee][1];
                    if (ntri[na] == vb and ntri[nb] == va) {
                        ne_geom = ee;
                        break;
                    }
                }
                std.debug.assert(ne_geom != 3);
                std.debug.assert(ne_geom == ne);
                // inverse perm check
                const pf = self.face_perm[f][e];
                const pb = self.face_perm[nf][ne];
                std.debug.assert(pb[pf[0]] == 0 and pb[pf[1]] == 1 and pb[pf[2]] == 2);
            }
        }
    }

    inline fn dirCubeStep(self: *const Grid, face: usize, d: u8) [3]isize {
        // Define step in LABEL basis then map to ACTUAL basis via A_f (label->actual)
        const Af = self.face_axes[face];
        var out: [3]isize = .{ 0, 0, 0 };
        out[Af[0]] = DIR_BARY_LABEL[d][0];
        out[Af[1]] = DIR_BARY_LABEL[d][1];
        out[Af[2]] = DIR_BARY_LABEL[d][2];
        return out;
    }

    inline fn pickViolatedEdgePost(self: *const Grid, face: usize, uvw1: [3]isize, _: [3]isize) ?u8 {
        // Low-bound only: choose seam only if a LABEL component is strictly negative.
        // Compute label-space components in A,B,C order deterministically.
        const Af = self.face_axes[face]; // label -> actual
        const uA = uvw1[Af[0]];
        const uB = uvw1[Af[1]];
        const uC = uvw1[Af[2]];
        if (uA < 0) return 0; // A<0 -> edge e=0 (B,C)
        if (uB < 0) return 1; // B<0 -> edge e=1 (C,A)
        if (uC < 0) return 2; // C<0 -> edge e=2 (A,B)
        // No seam to cross
        return null;
    }

    inline fn findDirByDelta(self: *const Grid, face: usize, delta_actual: [3]isize) ?u8 {
        // delta must be one of the 6 unit cube steps in actual basis (sum=0, components in {-1,0,1})
        std.debug.assert(delta_actual[0] + delta_actual[1] + delta_actual[2] == 0);
        var k: u8 = 0;
        while (k < 6) : (k += 1) {
            const cand = self.dirCubeStep(face, k);
            if (cand[0] == delta_actual[0] and cand[1] == delta_actual[1] and cand[2] == delta_actual[2]) {
                return k;
            }
        }
        return null;
    }

    pub inline fn isPentaAt(self: *const Grid, face: usize, q: isize, r: isize) bool {
        return self.isPentaFace(face, q, r);
    }

    pub inline fn isPentaIndex(self: *const Grid, idx: usize) bool {
        inline for (0..12) |gv| {
            if (self.penta_indices[gv] != 0 and self.penta_indices[gv] == idx) return true;
        }
        return false;
    }

    inline fn stepPermuteExact(self: *const Grid, face: usize, q: isize, r: isize, dir: u8) Coord {
        _ = .{ self, face, q, r, dir };
        legacyRemoved("stepPermuteExact removed; use prebaked TileId neighbors.");
        unreachable;
    }

    // normalize the "odd" edge coordinate for position extraction:
    // clamp to [N+1 .. 2N-1] and force odd parity deterministically.
    inline fn normalizeEdgeOdd(odd_in: isize, N: isize) isize {
        var odd = odd_in;
        const lo = N + 1;
        const hi = 2 * N - 1;
        if (odd < lo) odd = lo;
        if (odd > hi) odd = hi;
        if ((odd & 1) == 0) {
            // bias inward, preserving range. prefer predecessor when possible.
            if (odd > lo) {
                odd -= 1;
            } else {
                odd += 1; // lo was even; move up to first valid odd
            }
        }
        return odd;
    }

    inline fn edgeIndexForPair(a: usize, b: usize) usize {
        const x = if (a < b) a else b;
        const y = if (a < b) b else a;
        if (x == 0 and y == 1) return 0;
        if (x == 1 and y == 2) return 1;
        return 2; // {0,2}
    }

    // Oracle: return nearest global icosahedron vertex index (0..11) if within tolerance.
    pub inline fn getGlobalVertex(_g: *const Grid, face: usize, q: isize, r: isize) ?usize {
        // Use corrected spherical position
        const p = _g.faceCenterPosition(q, r, face);
        // Tolerance tuned for normalized unit sphere coordinates
        const tol: f64 = 1e-6;
        var best_idx: usize = 0;
        var best_d2: f64 = std.math.inf(f64);
        var vi: usize = 0;
        while (vi < icosa_vertices.len) : (vi += 1) {
            const v = icosa_vertices[vi];
            const dx = p.x - v.x;
            const dy = p.y - v.y;
            const dz = p.z - v.z;
            const d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best_d2) {
                best_d2 = d2;
                best_idx = vi;
            }
        }
        if (best_d2 <= tol) return best_idx;
        return null;
    }

    // Debug filters for isolating seam logs in tests; set to some value to enable
    const DBG_SEAM_FACE: ?usize = null;
    const DBG_SEAM_EI: ?usize = null;
    const DBG_SEAM_NF: ?usize = null;

    fn checkInboundUniqueness(self: *const Grid) bool {
        var ok = true;
        var nf: usize = 0;
        while (nf < 20) : (nf += 1) {
            var seen: [3]bool = .{ false, false, false };
            var f: usize = 0;
            while (f < 20) : (f += 1) {
                var e: usize = 0;
                while (e < 3) : (e += 1) {
                    if (self.face_neighbors[f][e][0] - 1 == nf) {
                        const ne = self.face_neighbors[f][e][1] - 1;
                        if (ne < 3) {
                            if (seen[ne]) ok = false;
                            seen[ne] = true;
                        } else {
                            ok = false;
                        }
                    }
                }
            }
            if (!(seen[0] and seen[1] and seen[2])) ok = false;
        }
        return ok;
    }

    fn countSetFaces(self: *const Grid) usize {
        const AX_UNSET: u8 = 255;
        var c: usize = 0;
        var i: usize = 0;
        while (i < 20) : (i += 1) {
            const P = self.face_axes[i];
            if (P[0] != AX_UNSET and P[1] != AX_UNSET and P[2] != AX_UNSET) c += 1;
        }
        return c;
    }

    inline fn expectedPnFromSeam(Pf: [3]u8, perm_fwd: [3]u8) [3]u8 {
        // Constraint: Pn[perm[j]] = Pf[j]
        var Pn: [3]u8 = .{ 0, 0, 0 };
        Pn[perm_fwd[0]] = Pf[0];
        Pn[perm_fwd[1]] = Pf[1];
        Pn[perm_fwd[2]] = Pf[2];
        return Pn;
    }

    inline fn expectedPermFromAxes(Pf: [3]u8, Pn: [3]u8) [3]u8 {
        var inv: [3]u8 = .{ 0, 0, 0 };
        inv[Pf[0]] = 0;
        inv[Pf[1]] = 1;
        inv[Pf[2]] = 2;
        var perm: [3]u8 = .{ 0, 0, 0 };
        perm[0] = inv[Pn[0]];
        perm[1] = inv[Pn[1]];
        perm[2] = inv[Pn[2]];
        return perm;
    }

    inline fn derivePermFromAxes(Pf: [3]u8, Pn_consensus: [3]u8) [3]u8 {
        return expectedPermFromAxes(Pf, Pn_consensus);
    }

    fn enqueueIfUnset(self: *Grid, q: *std.ArrayListUnmanaged(usize), f: usize) void {
        if (f >= 20) return;
        const P = self.face_axes[f];
        if (P[0] == 255 or P[1] == 255 or P[2] == 255) {
            _ = q.append(self.allocator, f) catch {};
        }
    }

    fn retargetInboundSeams(self: *Grid, nf: usize, old_ne: usize, new_src_f: usize, new_src_e: usize) void {
        if (nf >= 20 or old_ne >= 3) return;
        var g: usize = 0;
        while (g < 20) : (g += 1) {
            var ge: usize = 0;
            while (ge < 3) : (ge += 1) {
                const t_nf = self.face_neighbors[g][ge][0] - 1;
                const t_ne = self.face_neighbors[g][ge][1] - 1;
                if (t_nf == nf and t_ne == old_ne) {
                    if (!(g == new_src_f and ge == new_src_e)) {
                        // Clear stale mapping; keep a harmless placeholder
                        self.face_perm[g][ge] = .{ 0, 1, 2 };
                        self.face_neighbors[g][ge] = .{ nf + 1, old_ne + 1 };
                    }
                }
            }
        }
    }

    fn propagateAxesNoAssert(self: *Grid) void {
        // Seed Pf from face 0 and propagate without asserting on conflicts.
        const AX_UNSET: u8 = 255;
        var i: usize = 0;
        while (i < 20) : (i += 1) self.face_axes[i] = .{ AX_UNSET, AX_UNSET, AX_UNSET };
        self.face_axes[0] = .{ 0, 1, 2 };
        var q: [20]usize = undefined;
        var head: usize = 0;
        var tail: usize = 0;
        q[tail] = 0;
        tail += 1;
        while (head < tail) {
            const f = q[head];
            head += 1;
            const Pf = self.face_axes[f];
            if (Pf[0] == AX_UNSET) continue;
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf0 = self.face_neighbors[f][e][0];
                const ne0 = self.face_neighbors[f][e][1];
                if (nf0 == 0 or ne0 == 0) continue; // skip unset slots
                const nf = nf0 - 1;
                if (nf >= 20) continue;
                // Build expected neighbor axes purely by label-equality (standard basis propagation)
                const Vf = self.normalized_faces[f];
                const Vn = self.normalized_faces[nf];
                var expect: [3]u8 = .{ AX_UNSET, AX_UNSET, AX_UNSET };
                var used_mask: u8 = 0;
                var jj: usize = 0;
                while (jj < 3) : (jj += 1) {
                    // neighbor local vertex j should inherit axis from matching label on f
                    const label = Vn[jj];
                    var ii: usize = 0;
                    var found = false;
                    while (ii < 3) : (ii += 1) {
                        if (Vf[ii] == label) {
                            expect[jj] = Pf[ii];
                            const bit: u8 = @intCast(Pf[ii] & 3);
                            const mask_bit: u8 = if (bit == 0) 1 else if (bit == 1) 2 else 4;
                            used_mask |= mask_bit;
                            found = true;
                            break;
                        }
                    }
                    // It's okay if not found (the third vertex not shared); fill later
                }
                // Fill any unset slot with the remaining std axis
                const rem: u8 = if ((used_mask & 1) == 0) 0 else if ((used_mask & 2) == 0) 1 else 2;
                jj = 0;
                while (jj < 3) : (jj += 1) {
                    if (expect[jj] == AX_UNSET) {
                        expect[jj] = rem;
                        break;
                    }
                }
                const Pn = self.face_axes[nf];
                if (Pn[0] == AX_UNSET and Pn[1] == AX_UNSET and Pn[2] == AX_UNSET) {
                    self.face_axes[nf] = expect;
                    q[tail] = nf;
                    tail += 1;
                }
            }
        }
    }
    inline fn propagateAxesFromLabels(self: *Grid) void {
        // Alias to the label-only propagation with identity seed on face 0
        self.propagateAxesNoAssert();
    }

    inline fn toBaryFace(self: *const Grid, face: usize, q: isize, r: isize) [3]isize {
        const uvw = toBary(@intCast(self.size), q, r);
        const ax = self.face_axes[face];
        return .{ uvw[ax[0]], uvw[ax[1]], uvw[ax[2]] };
    }

    inline fn fromBaryFace(self: *const Grid, face: usize, uvw_face: [3]isize) QR {
        const ax = self.face_axes[face];
        var inv: [3]u8 = .{ 0, 0, 0 };
        inv[ax[0]] = 0;
        inv[ax[1]] = 1;
        inv[ax[2]] = 2;
        const uvw_std = .{ uvw_face[inv[0]], uvw_face[inv[1]], uvw_face[inv[2]] };
        return fromBary(@intCast(self.size), uvw_std);
    }

    inline fn isFaceHexValid(_: *const Grid, uvw_face: [3]isize, N: isize) bool {
        const maxB: isize = 2 * N;
        const eps: isize = 1; // tolerant by one integer to absorb rounding jitter at seams
        // box with tolerance
        if (uvw_face[0] < -eps or uvw_face[1] < -eps or uvw_face[2] < -eps) return false;
        if (uvw_face[0] > maxB + eps or uvw_face[1] > maxB + eps or uvw_face[2] > maxB + eps) return false;
        // triangle (u-w)<=N, (v-u)<=N, (w-v)<=N with tolerance
        if (uvw_face[0] - uvw_face[2] > N + eps) return false;
        if (uvw_face[1] - uvw_face[0] > N + eps) return false;
        if (uvw_face[2] - uvw_face[1] > N + eps) return false;
        return true;
    }

    inline fn takePlaneCandidate(num: isize, den: isize, kind: u8, idx: usize, best_num_ptr: *isize, best_den_ptr: *isize, any_ptr: *bool, best_kind_ptr: *u8, best_idx_ptr: *usize) void {
        if (den <= 0 or num < 0) return;
        if (!any_ptr.* or (num * best_den_ptr.*) < (best_num_ptr.* * den)) {
            any_ptr.* = true;
            best_kind_ptr.* = kind;
            best_idx_ptr.* = idx;
            best_num_ptr.* = num;
            best_den_ptr.* = den;
        }
    }

    inline fn clamp01(x: isize, lo: isize, hi: isize) isize {
        return if (x < lo) lo else if (x > hi) hi else x;
    }

    inline fn iabs(x: isize) isize {
        return if (x < 0) -x else x;
    }

    inline fn argmin3(a: isize, b: isize, c: isize) usize {
        if (a <= b and a <= c) return 0;
        if (b <= a and b <= c) return 1;
        return 2;
    }

    inline fn argmax3(a: isize, b: isize, c: isize) usize {
        if (a >= b and a >= c) return 0;
        if (b >= a and b >= c) return 1;
        return 2;
    }

    inline fn ceil_half_pos(e: isize) isize {
        return @divTrunc(e + 1, 2);
    }

    inline fn enforceHexWindow(uvw_in: [3]isize, N: isize) [3]isize {
        var u: isize = uvw_in[0];
        var v: isize = uvw_in[1];
        var w: isize = uvw_in[2];
        // Fix sum to exactly 3N by nudging smallest/largest as needed.
        var d: isize = 3 * N - (u + v + w);
        while (d != 0) {
            if (d > 0) {
                const t = argmin3(u, v, w);
                switch (t) {
                    0 => u += 1,
                    1 => v += 1,
                    else => w += 1,
                }
                d -= 1;
            } else {
                const t = argmax3(u, v, w);
                switch (t) {
                    0 => u -= 1,
                    1 => v -= 1,
                    else => w -= 1,
                }
                d += 1;
            }
        }
        // Symmetric deviation-space shrink: move extreme deviations toward center until spread <= N
        var d0: isize = u - N;
        var d1: isize = v - N;
        var d2: isize = w - N;
        // Invariants: d0 + d1 + d2 == 0 always holds
        while (true) {
            // find indices of min and max
            var lo_idx: usize = 0;
            var hi_idx: usize = 0;
            var lo_val: isize = d0;
            var hi_val: isize = d0;
            if (d1 < lo_val) {
                lo_val = d1;
                lo_idx = 1;
            }
            if (d2 < lo_val) {
                lo_val = d2;
                lo_idx = 2;
            }
            if (d1 > hi_val) {
                hi_val = d1;
                hi_idx = 1;
            }
            if (d2 > hi_val) {
                hi_val = d2;
                hi_idx = 2;
            }

            const spread: isize = hi_val - lo_val;
            if (spread <= N) break;
            const t: isize = ceil_half_pos(spread - N);
            switch (lo_idx) {
                0 => d0 += t,
                1 => d1 += t,
                else => d2 += t,
            }
            switch (hi_idx) {
                0 => d0 -= t,
                1 => d1 -= t,
                else => d2 -= t,
            }
        }
        u = N + d0;
        v = N + d1;
        w = N + d2;
        return .{ u, v, w };
    }

    // Helper function to hash coordinates for map lookup
    pub fn hashCoord(q: isize, r: isize, face: usize) u64 {
        const sz: i64 = 1 << 20;
        const q_u: u64 = @intCast(q + sz);
        const r_u: u64 = @intCast(r + sz);
        return (@as(u64, @intCast(face)) << 42) | (q_u << 21) | r_u;
    }

    pub fn unhashCoord(coord_hash: u64) struct { q: isize, r: isize, face: usize } {
        const face: usize = @intCast(coord_hash >> 42);
        const q_u: u64 = (coord_hash >> 21) & ((@as(u64, 1) << 21) - 1);
        const r_u: u64 = coord_hash & ((@as(u64, 1) << 21) - 1);
        const sz: i64 = 1 << 20;
        const q: isize = @intCast(@as(i64, @intCast(q_u)) - sz);
        const r: isize = @intCast(@as(i64, @intCast(r_u)) - sz);
        return .{ .q = q, .r = r, .face = face };
    }

    pub fn faceCenterPosition(self: *const Grid, q: isize, r: isize, face: usize) Vec3 {
        const N: isize = @intCast(self.size);
        const uvw_std = toBary(N, q, r);
        const ax = self.face_axes[face];
        const uvw_face = .{ uvw_std[ax[0]], uvw_std[ax[1]], uvw_std[ax[2]] };
        const tri = icosa_face_vertices[face]; // 1-based vertex ids
        const a = icosa_vertices[tri[0] - 1];
        const b = icosa_vertices[tri[1] - 1];
        const c = icosa_vertices[tri[2] - 1];
        const denom: f64 = @as(f64, @floatFromInt(@as(i64, @intCast(3 * N))));
        const wa: f64 = @as(f64, @floatFromInt(@as(i64, @intCast(uvw_face[0])))) / denom;
        const wb: f64 = @as(f64, @floatFromInt(@as(i64, @intCast(uvw_face[1])))) / denom;
        const wc: f64 = @as(f64, @floatFromInt(@as(i64, @intCast(uvw_face[2])))) / denom;
        var px = wa * a.x + wb * b.x + wc * c.x;
        var py = wa * a.y + wb * b.y + wc * c.y;
        var pz = wa * a.z + wb * b.z + wc * c.z;
        const len = std.math.sqrt(px * px + py * py + pz * pz);
        if (len != 0) {
            px /= len;
            py /= len;
            pz /= len;
        }
        return .{ .x = px, .y = py, .z = pz };
    }

    // Helper function to check if a hex coordinate is valid for a given size
    pub fn isValid(self: *const Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const N = @as(isize, @intCast(self.size));
        return iabs(q) <= N and iabs(r) <= N and iabs(s) <= N;
    }

    // Helper function to check if a tile is on the main icosahedron edge
    pub fn isEdge(self: *const Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const size_i = @as(isize, @intCast(self.size));
        return q - s == size_i or r - q == size_i or s - r == size_i;
    }

    // Helper function to check if a tile is a pentagon corner (face-local)
    pub fn isPenta(self: *const Grid, q: isize, r: isize) bool {
        // This function remains for existing call sites but is ambiguous without face.
        // Prefer using isPentaFace where possible.
        const face: usize = 0;
        return self.isPentaFace(face, q, r);
    }

    inline fn isPentaFace(self: *const Grid, face: usize, q: isize, r: isize) bool {
        const N: isize = @intCast(self.size);
        const uvw = self.toBaryFace(face, q, r);
        var zero_count: usize = 0;
        var has_2N = false;
        var has_N = false;
        var i: usize = 0;
        while (i < 3) : (i += 1) {
            if (uvw[i] == 0) zero_count += 1 else if (uvw[i] == 2 * N) has_2N = true else if (uvw[i] == N) has_N = true;
        }
        return zero_count == 1 and has_2N and has_N;
    }

    // Helper function to get which edge a coordinate is on (0..2)
    fn getEdgeIndex(self: *Grid, q: isize, r: isize) ?usize {
        const N: isize = @intCast(self.size);
        const s = -(q + r);
        if (q - s == N) return 0;
        if (r - q == N) return 1;
        if (s - r == N) return 2;
        return null;
    }

    // Test helpers (expose private edge helpers to tests)

    pub fn testEdgeAliasPair(self: *Grid, q: isize, r: isize, face: usize) ?struct { a: usize, b: usize } {
        if (self.isPenta(q, r)) return null;
        const ei = self.faceLocalEdgeIndex(face, q, r) orelse return null;
        const edge_len: usize = if (self.size > 0) self.size - 1 else 0;
        if (edge_len == 0) return null;
        const pos = self.getEdgePositionFace(q, r, face, ei);
        if (pos >= edge_len) return null;
        const idx_a = self.edge_indices[face][ei][pos];
        const map = self.neighborEdgeAndReversal(face, ei);
        const rev_pos = if (map.rev) (edge_len - 1) - pos else pos;
        if (map.nf >= 20 or map.ne >= 3 or rev_pos >= self.edge_indices[map.nf][map.ne].len) return null;
        const idx_b = self.edge_indices[map.nf][map.ne][rev_pos];
        return .{ .a = idx_a, .b = idx_b };
    }

    // Test helper: expose neighbor edge mapping for tests
    pub fn testNeighborEdgeMap(self: *const Grid, face: usize, ei: usize) EdgeMap {
        const map = self.neighborEdgeAndReversal(face, ei);
        return .{ .nf = map.nf, .ne = map.ne, .rev = map.rev };
    }

    // Test helper: expose enforceHexWindow for property-based tests
    pub fn testEnforceHexWindow(uvw_in: [3]isize, N: isize) [3]isize {
        return enforceHexWindow(uvw_in, N);
    }

    // Test helper: expose neighbor lookup (built at init)
    pub fn faceNeighborFace(self: *const Grid, face: usize, ei: usize) usize {
        return self.face_neighbors[face][ei][0] - 1;
    }

    // Test helper: expose exact stepper (preferred)

    // Test helper: expose face-aware edge position
    pub fn testGetEdgePositionFace(self: *const Grid, q: isize, r: isize, face: usize, edge: usize) usize {
        return self.getEdgePositionFace(q, r, face, edge);
    }

    // Test helper: expose face-aware barycentric conversion
    pub fn testToBaryFace(self: *const Grid, face: usize, q: isize, r: isize) [3]isize {
        return self.toBaryFace(face, q, r);
    }

    pub fn testFromBaryFace(self: *const Grid, face: usize, uvw_face: [3]isize) QR {
        return self.fromBaryFace(face, uvw_face);
    }
    // Test helper: map label-basis uvw to actual-basis uvw on a face
    pub fn testLabelToActual(self: *const Grid, face: usize, uvw_lbl: [3]isize) [3]isize {
        return applyPerm3(uvw_lbl, self.face_axes[face]);
    }
    // Test helper: expose seam-pick decision (optional)
    pub fn testPickViolatedEdge(self: *const Grid, face: usize, uvw1_actual: [3]isize) ?u8 {
        return self.pickViolatedEdgePost(face, uvw1_actual, .{ 0, 0, 0 });
    }
    // Test helper: exact stepper
    pub fn testExactStep(self: *const Grid, face: usize, q: isize, r: isize, dir: u8) Coord {
        _ = .{ self, face, q, r, dir };
        legacyRemoved("testExactStep removed; use prebaked TileId neighbors.");
        unreachable;
    }

    // Test helper: dir step delta in actual basis
    pub fn testDirCubeStep(self: *const Grid, face: usize, dir: u8) [3]isize {
        return self.dirCubeStep(face, dir);
    }
    // Test helper: find dir by exact delta
    pub fn testFindDirByDelta(self: *const Grid, face: usize, delta: [3]isize) ?u8 {
        return self.findDirByDelta(face, delta);
    }
    // Test helper: read a face's label at local index k
    pub fn testGetFaceLabel(self: *const Grid, face: usize, k: usize) u8 {
        return self.normalized_faces[face][k];
    }
    // Test helper: expose vertex owner for a global label
    pub fn testVertexOwner(self: *const Grid, label: u8) usize {
        return self.vertex_owner[label];
    }

    fn sortRing(self: *const Grid, center: usize, ids: []usize) void {
        const c = self.coords[center];
        const C = self.faceCenterPosition(c.q, c.r, c.face);
        const Cx = C.x;
        const Cy = C.y;
        const Cz = C.z;
        var t1x: f64 = -Cy;
        var t1y: f64 = Cx;
        var t1z: f64 = 0.0;
        var t1l = @sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
        if (t1l < 1e-12) {
            t1x = 0.0;
            t1y = -Cz;
            t1z = Cy;
            t1l = @sqrt(t1x * t1x + t1y * t1y + t1z * t1z);
        }
        const ux = t1x / t1l;
        const uy = t1y / t1l;
        const uz = t1z / t1l;
        const vx = Cy * uz - Cz * uy;
        const vy = Cz * ux - Cx * uz;
        const vz = Cx * uy - Cy * ux;
        // Compute angles for each id and selection-sort by angle (ring size <= 6)
        var angles: [6]f64 = undefined;
        var k: usize = 0;
        while (k < ids.len) : (k += 1) {
            const cc = self.coords[ids[k]];
            const P = self.faceCenterPosition(cc.q, cc.r, cc.face);
            const ax = P.x - Cx;
            const ay = P.y - Cy;
            const az = P.z - Cz;
            const axp = ax * ux + ay * uy + az * uz;
            const ayp = ax * vx + ay * vy + az * vz;
            angles[k] = std.math.atan2(ayp, axp);
        }
        var i: usize = 0;
        while (i + 1 < ids.len) : (i += 1) {
            var min_j = i;
            var j: usize = i + 1;
            while (j < ids.len) : (j += 1) {
                if (angles[j] < angles[min_j]) min_j = j;
            }
            if (min_j != i) {
                const tmp_id = ids[i];
                ids[i] = ids[min_j];
                ids[min_j] = tmp_id;
                const tmp_a = angles[i];
                angles[i] = angles[min_j];
                angles[min_j] = tmp_a;
            }
        }
    }

    // Position 0..(N-2) along interior of an edge (exclude pentagon corners)
    fn getEdgePosition(self: *Grid, q: isize, r: isize, edge: usize) usize {
        std.debug.assert(self.isEdge(q, r));
        std.debug.assert(!self.isPenta(q, r));
        const N: isize = @intCast(self.size);
        const uvw = toBary(N, q, r);
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        const b_local = pairs[edge][1];
        // derive from *current* edge bary; this is the true first-hit for edge tiles
        const odd = normalizeEdgeOdd(uvw[b_local], N);
        var pos = @divTrunc(odd - (N + 1), 2);
        if (pos < 0) pos = 0;
        const L: isize = N - 1;
        if (pos > L - 1) pos = L - 1;
        return @intCast(pos);
    }

    fn getEdgePositionFace(self: *const Grid, q: isize, r: isize, face: usize, edge: usize) usize {
        // tolerant: compute/clamp; caller should gate with faceLocalEdgeIndex
        const N: isize = @intCast(self.size);
        const uvw = self.toBaryFace(face, q, r);
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        const b_local = pairs[edge][1];
        const odd = normalizeEdgeOdd(uvw[b_local], N);
        var pos = @divTrunc(odd - (N + 1), 2);
        if (pos < 0) pos = 0;
        const L: isize = N - 1;
        if (pos > L - 1) pos = L - 1;
        return @intCast(pos);
    }

    fn neighborEdgeAndReversal(self: *const Grid, face: usize, ei: usize) EdgeMap {
        // neighbor face for this edge
        const nf = self.face_neighbors[face][ei][0] - 1;

        // Current face's shared edge vertices (local indices 0..2)
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        const a_local = pairs[ei][0];
        const b_local = pairs[ei][1];

        // Map them to neighbor-local indices via face_perm (current_idx -> neighbor_idx)
        const perm = self.face_perm[face][ei];
        const ia: usize = perm[a_local];
        const ib: usize = perm[b_local];

        // Neighbor edge index is the edge whose pair is {ia, ib}
        const ne: usize = edgeIndexForPair(ia, ib);

        // Define canonical forward as +1 mod 3
        const next = [_]usize{ 1, 2, 0 };
        const cur_forward = (b_local == next[a_local]);
        const nbr_forward = (ib == next[ia]);
        const rev = (cur_forward != nbr_forward);

        if (@import("builtin").is_test) {
            const match_face = if (DBG_SEAM_FACE) |fsel| (face == fsel) else false;
            const match_ei = if (DBG_SEAM_EI) |esel| (ei == esel) else true;
            const match_nf = if (DBG_SEAM_NF) |nsel| (nf == nsel) else true;
            if (match_face and match_ei and match_nf) {
                const Pf_dbg = self.face_axes[face];
                const Pn_dbg = self.face_axes[nf];
                std.debug.print(
                    "EDGE MAP: face={d} ei={d} -> nf={d} ia={d} ib={d} ne={d} rev={any} Pf={d},{d},{d} Pn={d},{d},{d}\n",
                    .{ face, ei, nf, ia, ib, ne, rev, Pf_dbg[0], Pf_dbg[1], Pf_dbg[2], Pn_dbg[0], Pn_dbg[1], Pn_dbg[2] },
                );
            }
        }
        return .{ .nf = nf, .ne = ne, .rev = rev };
    }

    inline fn faceLocalEdgeIndex(self: *const Grid, face: usize, q: isize, r: isize) ?usize {
        const uvw = self.toBaryFace(face, q, r);
        var zero: usize = 3;
        var i: usize = 0;
        while (i < 3) : (i += 1) {
            if (uvw[i] == 0) {
                zero = i;
                break;
            }
        }
        if (zero == 3) return null;
        return edgeIndexFromOppositeComponent(zero);
    }

    pub fn testFaceLocalEdgeIndex(self: *const Grid, face: usize, q: isize, r: isize) ?usize {
        return self.faceLocalEdgeIndex(face, q, r);
    }

    // Resolve index for pentagon tiles (12 special vertices of icosahedron)
    fn resolveIndexForPenta(self: *Grid, q: isize, r: isize, face: usize, current_index: usize) struct { usize, usize } {
        const N: isize = @intCast(self.size);
        if (face >= 20) return .{ current_index, current_index + 1 };
        // Determine local corner by which barycentric component equals 2N
        const uvw = self.toBaryFace(face, q, r);
        var local: usize = 0; // 0..2
        if (uvw[0] == 2 * N) local = 0 else if (uvw[1] == 2 * N) local = 1 else local = 2;
        const fv = icosa_face_vertices[face];
        const global_1 = fv[local];
        const gv = global_1 - 1; // 0..11
        if (gv >= self.penta_indices.len) return .{ current_index, current_index + 1 };
        if (self.penta_indices[gv] != UNSET) return .{ self.penta_indices[gv], current_index };
        self.penta_indices[gv] = current_index;
        return .{ current_index, current_index + 1 };
    }

    // Resolve index for edge tiles (tiles on the edges between faces)
    fn resolveIndexForEdge(self: *Grid, q: isize, r: isize, face: usize, current_index: usize) struct { usize, usize } {
        if (self.isPenta(q, r)) return .{ current_index, current_index + 1 };
        const ei = self.faceLocalEdgeIndex(face, q, r) orelse return .{ current_index, current_index + 1 };
        const edge_len: usize = if (self.size > 0) self.size - 1 else 0;
        if (edge_len == 0) return .{ current_index, current_index + 1 };
        const pos = self.getEdgePositionFace(q, r, face, ei);
        if (pos >= edge_len) return .{ current_index, current_index + 1 };
        if (self.edge_indices[face][ei][pos] != UNSET) return .{ self.edge_indices[face][ei][pos], current_index };
        self.edge_indices[face][ei][pos] = current_index;
        const map = self.neighborEdgeAndReversal(face, ei);
        const rev_pos = if (map.rev) (edge_len - 1) - pos else pos;
        if (debug_edge_logs < 30 and (face == 15 or face == 6)) {
            if (@import("builtin").is_test) {
                std.debug.print("edge alias: face={d} ei={d} pos={d} -> nf={d} ne={d} rev={any} rev_pos={d} q={d} r={d}\n", .{ face, ei, pos, map.nf, map.ne, map.rev, rev_pos, q, r });
                debug_edge_logs += 1;
            }
        }
        if (map.nf < 20 and map.ne < 3 and rev_pos < self.edge_indices[map.nf][map.ne].len) {
            self.edge_indices[map.nf][map.ne][rev_pos] = current_index;
        }
        return .{ current_index, current_index + 1 };
    }

    // New neighbor API: use helpers and an external coord_to_index map
    pub fn populateNeighbors(self: *Grid, coord_to_index: *const std.AutoHashMap(u64, usize)) !void {
        if (self.neighbors.len > 0) self.allocator.free(self.neighbors);
        self.neighbors = try self.allocator.alloc(NeighborSet, self.tile_count);

        // Alias sets no longer needed with global canonicalization

        // Normalize homes from canonical map; no face-0 preference required

        // Precompute pentagon classification per tile index from penta_indices
        var is_penta_index = try self.allocator.alloc(bool, self.tile_count);
        defer self.allocator.free(is_penta_index);
        @memset(is_penta_index, false);
        {
            var gv: usize = 0;
            while (gv < self.penta_indices.len) : (gv += 1) {
                const idx = self.penta_indices[gv];
                if (idx != 0 and idx < self.tile_count) is_penta_index[idx] = true;
            }
            if (@import("builtin").is_test) {
                var pcount: usize = 0;
                var tgv: usize = 0;
                while (tgv < self.penta_indices.len) : (tgv += 1) {
                    if (self.penta_indices[tgv] != 0) pcount += 1;
                }
                // Standard icosa grids should have 12 pentagon indices
                std.debug.assert(pcount == 12);
            }
        }

        var i: usize = 0;
        while (i < self.tile_count) : (i += 1) {
            // initialize ring to self and reset count to 0 before filling
            @memset(&self.neighbors[i].ids, i);
            self.neighbors[i].count = 0;
            const cap: u8 = if (self.isPentaAt(self.coords[i].face, self.coords[i].q, self.coords[i].r)) 5 else 6;
            const home = self.coords[i];
            // Build ring from home alias using global canonicalization of neighbors
            var out_ids: [6]usize = undefined;
            var out_count: u8 = 0;
            var dir: u8 = 0;
            while (dir < 6 and out_count < cap) : (dir += 1) {
                const fwd_raw = self.stepPermuteExact(home.face, home.q, home.r, dir);
                const cc = self.canonicalCoord(fwd_raw.face, fwd_raw.q, fwd_raw.r);
                const h = hashCoord(cc.q, cc.r, cc.face);
                if (coord_to_index.get(h)) |j| {
                    if (j == i) continue;
                    var seen = false;
                    var s: usize = 0;
                    while (s < out_count) : (s += 1) {
                        if (out_ids[s] == j) {
                            seen = true;
                            break;
                        }
                    }
                    if (!seen) {
                        out_ids[out_count] = j;
                        out_count += 1;
                    }
                }
            }
            // Sort ring by angle around center and deduplicate
            self.sortRing(i, out_ids[0..out_count]);
            out_count = @intCast(dedupSorted(out_ids[0..out_count]));
            if (out_count > cap) out_count = cap;
            var w: usize = 0;
            while (w < out_count) : (w += 1) self.neighbors[i].ids[w] = out_ids[w];
            // Optional: pad remaining slots with self for fixed-size storage
            while (w < cap) : (w += 1) self.neighbors[i].ids[w] = i;
            // Store the actual neighbor count (not capacity)
            self.neighbors[i].count = @intCast(out_count);

            // One-time probe: log first non-penta ring with degree < 6 to isolate seam issues
            if (@import("builtin").is_test and TEST_VERBOSE and !is_penta_index[i] and self.neighbors[i].count < 6 and !test_probe_done) {
                test_probe_done = true;
                const center = home;
                std.debug.print("RING DROP: center i={d} (f={d},q={d},r={d}) count={d}\n", .{ i, center.face, center.q, center.r, self.neighbors[i].count });
                self.auditRingWithMap(coord_to_index, i);
                var dd: u8 = 0;
                while (dd < 6) : (dd += 1) {
                    self.auditSeamTransport(i, dd);
                }
            }
        }

        // No symmetry backfill: correctness must come from exact stepper + canonicalization
    }

    // Debug-only: ring audit using exact stepper and canonicalization
    fn auditRingWithMap(self: *Grid, coord_to_index: *const std.AutoHashMap(u64, usize), i: usize) void {
        if (!@import("builtin").is_test or !TEST_VERBOSE) return;
        const home = self.coords[i];
        var seen = std.AutoHashMap(usize, u8).init(self.allocator);
        defer seen.deinit();
        var d: u8 = 0;
        while (d < 6) : (d += 1) {
            const nb = self.stepPermuteExact(home.face, home.q, home.r, d);
            const can = self.canonicalCoord(nb.face, nb.q, nb.r);
            const h = hashCoord(can.q, can.r, can.face);
            const j_opt = coord_to_index.get(h);
            if (j_opt == null) {
                std.debug.print("RING AUDIT i={d} d={d}: unresolved -> (f={d},q={d},r={d})\n", .{ i, d, can.face, can.q, can.r });
                continue;
            }
            const j = j_opt.?;
            if (seen.contains(j)) {
                const prev_d = seen.get(j).?;
                std.debug.print("RING AUDIT i={d} DUP: d={d} and d={d} -> j={d} (f={d},q={d},r={d})\n", .{ i, d, prev_d, j, can.face, can.q, can.r });
            } else {
                _ = seen.put(j, d) catch {};
                std.debug.print("RING AUDIT i={d} d={d} -> j={d} (f={d},q={d},r={d})\n", .{ i, d, j, can.face, can.q, can.r });
            }
        }
    }

    // Debug-only: audit seam transport for a given tile/direction using exact integer logic
    fn auditSeamTransport(self: *Grid, i: usize, d: u8) void {
        if (!@import("builtin").is_test or !TEST_VERBOSE) return;
        const N: isize = @intCast(self.size);
        const R = 2 * N; // component bounds in uvw
        const h = self.coords[i];
        const uvw = self.toBaryFace(h.face, h.q, h.r);
        const duvw = self.dirCubeStep(h.face, d);
        const uvw1 = .{ uvw[0] + duvw[0], uvw[1] + duvw[1], uvw[2] + duvw[2] };
        const cn_opt = self.pickViolatedEdgePost(h.face, uvw1, duvw);
        if (cn_opt == null) return;
        const e = Grid.edgeIndexFromOppositeComponent(cn_opt.?);
        const nf = self.face_neighbors[h.face][e][0] - 1;
        const perm = self.face_perm[h.face][e];
        const uvw_nf = applyPerm3(uvw, perm);
        const duvw_nf = applyPerm3(duvw, perm);
        const uvw2 = .{ uvw_nf[0] + duvw_nf[0], uvw_nf[1] + duvw_nf[1], uvw_nf[2] + duvw_nf[2] };
        std.debug.assert(uvw2[0] >= 0 and uvw2[1] >= 0 and uvw2[2] >= 0 and uvw2[0] <= R and uvw2[1] <= R and uvw2[2] <= R);
        const nb_qr = self.fromBaryFace(nf, uvw2);
        // Find reverse dir on nf by matching -duvw_nf
        const rdu = -duvw_nf[0];
        const rdv = -duvw_nf[1];
        const rdw = -duvw_nf[2];
        var rev_d: u8 = 255;
        var k: u8 = 0;
        while (k < 6) : (k += 1) {
            const cd = self.dirCubeStep(nf, k);
            if (cd[0] == rdu and cd[1] == rdv and cd[2] == rdw) {
                rev_d = k;
                break;
            }
        }
        if (rev_d == 255) {
            std.debug.print("SEAM AUDIT i={d} d={d}: no reverse dir on nf={d} e={d} perm=({d},{d},{d})\n", .{ i, d, nf, e, perm[0], perm[1], perm[2] });
            return;
        }
        const back = self.stepPermuteExact(nf, nb_qr.q, nb_qr.r, rev_d);
        const back_can = self.canonicalCoord(back.face, back.q, back.r);
        const home_can = self.canonicalCoord(h.face, h.q, h.r);
        if (!self.coordEq(back_can, home_can)) {
            std.debug.print("SEAM AUDIT i={d} d={d}: nf={d} e={d} perm=({d},{d},{d}) nb=(f={d},q={d},r={d}) rev_d={d} -> back=(f={d},q={d},r={d}) (FAIL)\n", .{ i, d, nf, e, perm[0], perm[1], perm[2], nf, nb_qr.q, nb_qr.r, rev_d, back_can.face, back_can.q, back_can.r });
        } else {
            std.debug.print("SEAM AUDIT i={d} d={d}: nf={d} e={d} perm=({d},{d},{d}) rev_d={d} (OK)\n", .{ i, d, nf, e, perm[0], perm[1], perm[2], rev_d });
        }
    }

    const TEST_VERBOSE: bool = true;

    fn stepAcrossOrIn(self: *const Grid, q: isize, r: isize, face: usize, dir_idx: u8, _: *const anyopaque) Coord {
        return self.stepPermuteExact(face, q, r, dir_idx);
    }

    pub fn populateIndices(self: *Grid) !std.AutoHashMap(u64, usize) {
        _ = self;
        legacyRemoved("populateIndices removed; prebaked builder enumerates analytically.");
        unreachable;
    }

    pub fn generateMesh(self: *Grid) !void {
        if (self.vertices) |v| self.allocator.free(v);
        if (self.elevations) |e| self.allocator.free(e);
        if (self.indices) |i| self.allocator.free(i);
        self.vertices = null;
        self.elevations = null;
        self.indices = null;

        const n = self.tile_count;
        if (n == 0) return;

        const alloc = self.allocator;

        var verts = try alloc.alloc(f32, n * 3);
        errdefer alloc.free(verts);
        var elevs = try alloc.alloc(f32, n);
        errdefer alloc.free(elevs);

        var i: usize = 0;
        while (i < n) : (i += 1) {
            const c = self.coords[i];
            const pos = self.faceCenterPosition(c.q, c.r, c.face);
            const base = i * 3;
            verts[base + 0] = @floatCast(pos.x);
            verts[base + 1] = @floatCast(pos.y);
            verts[base + 2] = @floatCast(pos.z);
            elevs[i] = 0.0;
        }

        var idx_list: std.ArrayListUnmanaged(u32) = .{};
        errdefer idx_list.deinit(alloc);

        var a: usize = 0;
        while (a < n) : (a += 1) {
            const ring = self.neighbors[a];
            const m: u8 = ring.count;
            if (m < 2) continue;
            var j: u8 = 0;
            while (j < m) : (j += 1) {
                const b_idx = ring.ids[j];
                const c_idx = ring.ids[@intCast((@as(u16, j) + 1) % @as(u16, m))];
                if (b_idx == a or c_idx == a or b_idx == c_idx) continue;
                if (a < b_idx and a < c_idx) {
                    try idx_list.append(alloc, @intCast(a));
                    try idx_list.append(alloc, @intCast(b_idx));
                    try idx_list.append(alloc, @intCast(c_idx));
                }
            }
        }

        self.vertices = verts;
        self.elevations = elevs;
        self.indices = try idx_list.toOwnedSlice(alloc);
    }

    // Test helper: print full violation vector and seam decision for diagnostics
    pub fn testSeamDebug(self: *const Grid, q: isize, r: isize, face: usize, dir_idx: u8) void {
        _ = self.stepPermuteExact(face, q, r, dir_idx);
    }

    // Test helper: expose geometry-derived perm for seam (f,e)
    pub fn testBuildPermForEdge(self: *const Grid, f: usize, e: usize) [3]u8 {
        return self.buildPermForEdge(f, e);
    }

    pub fn testExpectedPermFromAxes(_: *const Grid, Pf: [3]u8, Pn: [3]u8) [3]u8 {
        return expectedPermFromAxes(Pf, Pn);
    }

    // Test helper: return a step result (q,r,face) triple
    pub fn testStep(self: *const Grid, q: isize, r: isize, face: usize, dir_idx: u8) struct { q: isize, r: isize, face: usize } {
        const c = self.stepPermuteExact(face, q, r, dir_idx);
        return .{ .q = c.q, .r = c.r, .face = c.face };
    }

    // Post-BFS alias completion for pentagon vertices: ensure all 5 face-aliases exist per global vertex

    // Debug: audit a single corner-turn table entry against topology
    pub fn debugCornerTurnEntry(self: *const Grid, f: usize, v_local: usize, eidx: usize) void {
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        const tri_f = icosa_face_vertices[f];
        const nf = self.face_neighbors[f][eidx][0] - 1;
        const ne = self.face_neighbors[f][eidx][1] - 1;
        const tri_nf = icosa_face_vertices[nf];
        const gv = tri_f[v_local];
        var v_nf: usize = 3;
        var i: usize = 0;
        while (i < 3) : (i += 1) {
            if (tri_nf[i] == gv) {
                v_nf = i;
                break;
            }
        }
        var e2_expected: usize = 3;
        var k: usize = 0;
        while (k < 3) : (k += 1) {
            const inc = (pairs[k][0] == v_nf) or (pairs[k][1] == v_nf);
            if (!inc or k == ne) continue;
            e2_expected = k;
            break;
        }
        const nf2_expected = if (e2_expected < 3) self.face_neighbors[nf][e2_expected][0] - 1 else 999;
        if (@import("builtin").is_test and TEST_VERBOSE) {
            std.debug.print("CTDBG f={d} v={d} e={d} -> nf={d} ne={d} | e2(exp={d}) nf2(exp={d})\n", .{
                f, v_local, eidx, nf, ne, e2_expected, nf2_expected,
            });
        }
    }

    // Helper: find face-local vertex index for a hex corner tile (argmax component in face bary)
    pub fn faceLocalVertexForCorner(self: *const Grid, f: usize, q: isize, r: isize) usize {
        const uvw = self.toBaryFace(f, q, r);
        var vmax: isize = -1;
        var vi: usize = 0;
        inline for (0..3) |i| {
            if (uvw[i] > vmax) {
                vmax = uvw[i];
                vi = i;
            }
        }
        return vi;
    }

    // Small helper: safely read neighbor (0-based) or return null if unset.
    inline fn getNeighbor0Based(self: *const Grid, face: usize, e: usize) ?struct { nf: usize, ne: usize } {
        const nei = self.face_neighbors[face][e];
        if (nei[0] == 0 or nei[1] == 0) return null;
        return .{ .nf = nei[0] - 1, .ne = nei[1] - 1 };
    }

    // Small helper: check whether axes are fully set for a face (no 255 sentinels).
    inline fn axesSet(self: *const Grid, face: usize) bool {
        const ax = self.face_axes[face];
        return (ax[0] < 3 and ax[1] < 3 and ax[2] < 3);
    }
};

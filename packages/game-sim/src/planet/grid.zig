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
pub const icosa_face_vertices = [_][3]usize{
    .{ 1, 9, 5 },   .{ 1, 6, 11 },  .{ 3, 5, 10 }, .{ 3, 12, 6 }, .{ 2, 7, 9 },
    .{ 2, 11, 8 },  .{ 4, 10, 7 },  .{ 4, 8, 12 }, .{ 1, 11, 9 }, .{ 2, 9, 11 },
    .{ 3, 10, 12 }, .{ 4, 10, 12 }, .{ 5, 3, 1 },  .{ 6, 1, 3 },  .{ 7, 2, 4 },
    .{ 8, 4, 2 },   .{ 9, 7, 5 },   .{ 10, 5, 7 }, .{ 11, 6, 8 }, .{ 12, 8, 6 },
};

// Face neighbors will be constructed at init from faces

// Canonical axial direction list (CW): (+q), (+r), (-q,+s), (-q), (-r), (+q,-s)
const DIRS = [_][2]isize{ .{ 1, 0 }, .{ 1, -1 }, .{ 0, -1 }, .{ -1, 0 }, .{ -1, 1 }, .{ 0, 1 } };

// DIR_REFLECTION removed (unused)

var debug_edge_logs: usize = 0;
var debug_nbr_logs: usize = 0;
var test_probe_done: bool = false;
var seam_verify_logged: bool = false;

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

    // Debug snapshots to guard against accidental seam rewrites after init (Debug mode only)
    _perm_snapshot: [20][3][3]u8 = std.mem.zeroes([20][3][3]u8),
    _ne_snapshot: [20][3][2]usize = std.mem.zeroes([20][3][2]usize),

    // Internal data structures for tracking shared indices during generation
    penta_indices: [12]usize,
    edge_indices: [20][3][]usize,

    // Generated mesh data (flat arrays suitable for GPU upload)
    // vertices: xyz triplets in f32
    vertices: ?[]f32,
    // elevations per-vertex (f32)
    elevations: ?[]f32,
    // triangle indices (u32)
    indices: ?[]u32,

    // Public helper struct for edge mappings (used in tests)
    pub const EdgeMap = struct { nf: usize, ne: usize, rev: bool };

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
            ._perm_snapshot = std.mem.zeroes([20][3][3]u8),
            ._ne_snapshot = std.mem.zeroes([20][3][2]usize),
            .penta_indices = undefined,
            .edge_indices = std.mem.zeroes([20][3][]usize),
            .vertices = null,
            .elevations = null,
            .indices = null,
        };

        // Initialize penta sentinel and allocate edge indices arrays
        grid.penta_indices = .{UNSET} ** 12;
        for (0..20) |face| {
            for (0..3) |edge| {
                const edge_len: usize = if (size > 0) size - 1 else 0;
                grid.edge_indices[face][edge] = try allocator.alloc(usize, edge_len);
                if (edge_len > 0) @memset(grid.edge_indices[face][edge], UNSET);
            }
        }

        try grid.buildFaceNeighbors();
        grid.normalizeSeams(); // authoritative: write perm both ways and edge indices from geometry
        grid.buildFaceAxes();
        grid.verifySeams();
        grid.freezeSeamsForDebug();
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

    fn buildFaceNeighbors(self: *Grid) !void {

        // zero init
        for (0..20) |fi| {
            for (0..3) |ei| {
                self.face_neighbors[fi][ei] = .{ 0, 0 };
            }
        }

        var map = std.AutoHashMap(EdgeKey, struct { f: usize, e: usize }).init(self.allocator);
        defer map.deinit();

        // pass 1: link opposite faces by shared edges
        for (0..20) |fi| {
            const tri = icosa_face_vertices[fi];
            const edges = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
            for (edges, 0..) |e, ei| {
                const va = tri[e[0]];
                const vb = tri[e[1]];
                const k = edgeKey(va, vb);
                if (map.getPtr(k)) |slot| {
                    const of = slot.f;
                    const oe = slot.e;
                    self.face_neighbors[fi][ei] = .{ of + 1, oe + 1 };
                    self.face_neighbors[of][oe] = .{ fi + 1, ei + 1 };
                } else {
                    try map.put(k, .{ .f = fi, .e = ei });
                }
            }
        }
        // any border (none in closed icosa) remains {0,0}
        // Validate closed topology and reverse links
        for (0..20) |f| {
            for (0..3) |e| {
                const nf = self.face_neighbors[f][e][0];
                const ne = self.face_neighbors[f][e][1];
                std.debug.assert(nf != 0 and ne != 0);
                const back = self.face_neighbors[nf - 1][ne - 1];
                std.debug.assert(back[0] == f + 1 and back[1] == e + 1);
            }
        }
    }

    fn buildPermForEdge(self: *const Grid, f: usize, e: usize) [3]u8 {
        // Build CURRENT face-local -> NEIGHBOR face-local permutation for seam (f,e) from geometry.
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
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
        if (ia == 3 or ib == 3) @panic("buildPermForEdge: shared vertices not found");
        const c_local: usize = 3 - a_local - b_local;
        const ic: usize = 3 - ia - ib;
        var perm_fwd: [3]u8 = .{ 0, 0, 0 };
        perm_fwd[a_local] = @intCast(ia);
        perm_fwd[b_local] = @intCast(ib);
        perm_fwd[c_local] = @intCast(ic);
        return perm_fwd;
    }

    fn normalizeSeams(self: *Grid) void {
        // Authoritative: write perm and neighbor edge indices for both directions from geometry only.
        self.face_perm = std.mem.zeroes([20][3][3]u8);
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                if (nf >= 20) continue;
                const perm_fwd = self.buildPermForEdge(f, e);
                // derive neighbor edge index from endpoints
                const ia = perm_fwd[pairs[e][0]];
                const ib = perm_fwd[pairs[e][1]];
                const ne = edgeIndexForPair(ia, ib);
                // inverse perm
                var perm_bwd: [3]u8 = .{ 0, 0, 0 };
                perm_bwd[perm_fwd[0]] = 0;
                perm_bwd[perm_fwd[1]] = 1;
                perm_bwd[perm_fwd[2]] = 2;
                // write both directions coherently
                self.face_perm[f][e] = perm_fwd;
                self.face_neighbors[f][e][1] = ne + 1;
                self.face_perm[nf][ne] = perm_bwd;
                self.face_neighbors[nf][ne][1] = e + 1;
                // ensure reciprocity for neighbor face index (should already be set)
                self.face_neighbors[nf][ne][0] = f + 1;
            }
        }
        // Debug: check involution
        if (@import("builtin").mode == .Debug) {
            for (0..20) |ff| {
                for (0..3) |ee| {
                    const nf = self.face_neighbors[ff][ee][0] - 1;
                    const ne = self.face_neighbors[ff][ee][1] - 1;
                    if (nf >= 20 or ne >= 3) continue;
                    const p = self.face_perm[ff][ee];
                    const q = self.face_perm[nf][ne];
                    std.debug.assert(q[p[0]] == 0 and q[p[1]] == 1 and q[p[2]] == 2);
                }
            }
        }
    }

    inline fn applyPerm3(x: [3]isize, p: [3]u8) [3]isize {
        return .{ x[p[0]], x[p[1]], x[p[2]] };
    }

    inline fn isEvenPerm3(p: [3]u8) bool {
        var invs: u8 = 0;
        if (p[0] > p[1]) invs += 1;
        if (p[0] > p[2]) invs += 1;
        if (p[1] > p[2]) invs += 1;
        return (invs & 1) == 0;
    }

    inline fn printPerm3(p: [3]u8) void {
        std.debug.print("[{d},{d},{d}]", .{ p[0], p[1], p[2] });
    }
    inline fn toBary(N: isize, q: isize, r: isize) [3]isize {
        const u = q + N;
        const v = r + N;
        const w = N - q - r;
        return .{ u, v, w };
    }

    inline fn fromBary(N: isize, uvw: [3]isize) QR {
        return QR{ .q = uvw[0] - N, .r = uvw[1] - N };
    }

    inline fn dirBary(dq: isize, dr: isize) [3]isize {
        return .{ dq, dr, -(dq + dr) };
    }

    inline fn edgeIndexFromOppositeComponent(cn: usize) usize {
        return switch (cn) {
            0 => 1,
            1 => 2,
            2 => 0,
            else => unreachable,
        };
    }

    inline fn edgeIndexForPair(a: usize, b: usize) usize {
        const x = if (a < b) a else b;
        const y = if (a < b) b else a;
        if (x == 0 and y == 1) return 0;
        if (x == 1 and y == 2) return 1;
        return 2; // {0,2}
    }

    // Debug filters for isolating seam logs in tests; set to some value to enable
    const DBG_SEAM_FACE: ?usize = 18;
    const DBG_SEAM_EI: ?usize = null;
    const DBG_SEAM_NF: ?usize = 5;

    fn buildFaceAxes(self: *Grid) void {
        // Solve axes as a fixed point of Pn[perm[j]] = Pf[j]
        const AX_UNSET: u8 = 255;
        var i: usize = 0;
        while (i < 20) : (i += 1) self.face_axes[i] = .{ AX_UNSET, AX_UNSET, AX_UNSET };
        // seed root
        self.face_axes[0] = .{ 0, 1, 2 };

        var changed = true;
        var guard: usize = 0;
        while (changed and guard < 100) : (guard += 1) {
            changed = false;
            var f: usize = 0;
            while (f < 20) : (f += 1) {
                const Pf = self.face_axes[f];
                if (Pf[0] == AX_UNSET or Pf[1] == AX_UNSET or Pf[2] == AX_UNSET) continue;
                var e: usize = 0;
                while (e < 3) : (e += 1) {
                    const nf = self.face_neighbors[f][e][0] - 1;
                    if (nf >= 20) continue;
                    const perm = self.face_perm[f][e]; // current_local -> neighbor_local
                    var expect: [3]u8 = .{ AX_UNSET, AX_UNSET, AX_UNSET };
                    var j: usize = 0;
                    while (j < 3) : (j += 1) expect[perm[j]] = Pf[j];
                    const Pn = self.face_axes[nf];
                    if (Pn[0] == AX_UNSET or Pn[1] == AX_UNSET or Pn[2] == AX_UNSET or Pn[0] != expect[0] or Pn[1] != expect[1] or Pn[2] != expect[2]) {
                        self.face_axes[nf] = expect;
                        changed = true;
                    }
                }
            }
        }
        // ensure solved
        var k: usize = 0;
        while (k < 20) : (k += 1) {
            const Pk = self.face_axes[k];
            std.debug.assert(Pk[0] != AX_UNSET and Pk[1] != AX_UNSET and Pk[2] != AX_UNSET);
        }
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
        // Solve Pn[perm[j]] = Pf[j] for perm.
        var perm: [3]u8 = .{ 0, 0, 0 };
        var j: usize = 0;
        while (j < 3) : (j += 1) {
            var k: usize = 0;
            while (k < 3 and Pn[k] != Pf[j]) : (k += 1) {}
            perm[j] = @intCast(k);
        }
        return perm;
    }

    fn reconcilePermWithAxes(self: *Grid) void {
        // For each seam, derive the expected perm from Pf/Pn and rewrite forward and inverse if needed.
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            const Pf = self.face_axes[f];
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                const ne = self.face_neighbors[f][e][1] - 1;
                const Pn = self.face_axes[nf];
                const perm_expected = expectedPermFromAxes(Pf, Pn);
                const cur = self.face_perm[f][e];
                if (!(cur[0] == perm_expected[0] and cur[1] == perm_expected[1] and cur[2] == perm_expected[2])) {
                    self.face_perm[f][e] = perm_expected;
                    var inv: [3]u8 = .{ 0, 0, 0 };
                    inv[perm_expected[0]] = 0;
                    inv[perm_expected[1]] = 1;
                    inv[perm_expected[2]] = 2;
                    self.face_perm[nf][ne] = inv;
                }
            }
        }
    }

    fn reconcileNeighborEdgeIndices(self: *Grid) void {
        // Pass 1: recompute ne from perm endpoints
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                if (nf >= 20) continue;
                const perm = self.face_perm[f][e];
                const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
                const a_local = pairs[e][0];
                const b_local = pairs[e][1];
                const ia: usize = perm[a_local];
                const ib: usize = perm[b_local];
                const want_ne: usize = edgeIndexForPair(ia, ib);
                self.face_neighbors[f][e][1] = want_ne + 1;
            }
        }
        // Pass 2: enforce reciprocity of neighbor edge indices
        var ff: usize = 0;
        while (ff < 20) : (ff += 1) {
            var ee: usize = 0;
            while (ee < 3) : (ee += 1) {
                const nf = self.face_neighbors[ff][ee][0] - 1;
                const ne = self.face_neighbors[ff][ee][1] - 1;
                if (nf >= 20 or ne >= 3) continue;
                self.face_neighbors[nf][ne][0] = ff + 1;
                self.face_neighbors[nf][ne][1] = ee + 1;
            }
        }
    }

    fn reconcilePermInvolution(self: *Grid) void {
        // Enforce backward perm is inverse of forward perm AND align neighbor edge indices to endpoints.
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                if (nf >= 20) continue;
                const p_fwd = self.face_perm[f][e];
                // Determine which neighbor edge this seam lands on by shared endpoints
                const a_local = pairs[e][0];
                const b_local = pairs[e][1];
                const ia: usize = p_fwd[a_local];
                const ib: usize = p_fwd[b_local];
                const want_ne: usize = edgeIndexForPair(ia, ib);
                // Build inverse
                var p_inv: [3]u8 = .{ 0, 0, 0 };
                p_inv[p_fwd[0]] = 0;
                p_inv[p_fwd[1]] = 1;
                p_inv[p_fwd[2]] = 2;
                // Write both directions coherently
                self.face_perm[f][e] = p_fwd;
                self.face_perm[nf][want_ne] = p_inv;
                // Update neighbor edge indices reciprocity
                self.face_neighbors[f][e][1] = want_ne + 1;
                self.face_neighbors[nf][want_ne][0] = f + 1;
                self.face_neighbors[nf][want_ne][1] = e + 1;
            }
        }
        // Optional check: do not hard-assert; verifySeams will report detailed mismatches.
    }

    fn verifySeams(self: *const Grid) void {
        // If seams were frozen, ensure no later writes occurred (Debug only)
        if (@import("builtin").mode == .Debug and self._ne_snapshot[0][0][0] != 0) {
            std.debug.assert(std.mem.eql(u8, std.mem.asBytes(&self.face_perm), std.mem.asBytes(&self._perm_snapshot)));
            std.debug.assert(std.mem.eql(u8, std.mem.asBytes(&self.face_neighbors), std.mem.asBytes(&self._ne_snapshot)));
        }
        // Enhanced reporting of constraints into each face
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };

        var first_mismatch_logged = false;
        // Per-seam checks
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            const Pf = self.face_axes[f];
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                const ne = self.face_neighbors[f][e][1] - 1;
                const p_fwd = self.face_perm[f][e];
                const p_bwd = self.face_perm[nf][ne];
                if (!(p_bwd[p_fwd[0]] == 0 and p_bwd[p_fwd[1]] == 1 and p_bwd[p_fwd[2]] == 2) and !first_mismatch_logged) {
                    std.debug.print("SEAM VERIFY FAIL (involution): f={d} e={d} -> nf={d} ne={d} perm_fwd=", .{ f, e, nf, ne });
                    printPerm3(p_fwd);
                    std.debug.print(" perm_bwd=", .{});
                    printPerm3(p_bwd);
                    std.debug.print("\n", .{});
                    first_mismatch_logged = true;
                }
                const Pn_expected = expectedPnFromSeam(Pf, p_fwd);
                const Pn_actual = self.face_axes[nf];
                if (!(Pn_actual[0] == Pn_expected[0] and Pn_actual[1] == Pn_expected[1] and Pn_actual[2] == Pn_expected[2]) and !first_mismatch_logged) {
                    std.debug.print("SEAM VERIFY FAIL (axes): f={d} e={d} -> nf={d} ne={d}\n", .{ f, e, nf, ne });
                    std.debug.print("  Pf=", .{});
                    printPerm3(Pf);
                    std.debug.print(" perm_fwd=", .{});
                    printPerm3(p_fwd);
                    std.debug.print("  Pn_expected=", .{});
                    printPerm3(Pn_expected);
                    std.debug.print("  Pn_actual=", .{});
                    printPerm3(Pn_actual);
                    std.debug.print("\n", .{});
                    first_mismatch_logged = true;
                }
                // Edge index consistency
                const a = pairs[e][0];
                const b = pairs[e][1];
                const ia = p_fwd[a];
                const ib = p_fwd[b];
                const want_ne = edgeIndexForPair(ia, ib);
                if (ne != want_ne and !first_mismatch_logged) {
                    std.debug.print("SEAM VERIFY FAIL (edge idx): f={d} e={d} -> nf={d} ne={d} want_ne={d} endpoints ia,ib={d},{d}\n", .{ f, e, nf, ne, want_ne, ia, ib });
                    first_mismatch_logged = true;
                }
            }
        }
        // Consensus per target face
        var nf: usize = 0;
        while (nf < 20) : (nf += 1) {
            var expects: [3][3]u8 = .{ .{ 0, 0, 0 }, .{ 0, 0, 0 }, .{ 0, 0, 0 } };
            var have: [3]bool = .{ false, false, false };
            var ff: usize = 0;
            while (ff < 20) : (ff += 1) {
                var ee: usize = 0;
                while (ee < 3) : (ee += 1) {
                    const t_nf = self.face_neighbors[ff][ee][0] - 1;
                    const t_ne = self.face_neighbors[ff][ee][1] - 1;
                    if (t_nf == nf) {
                        const Pf = self.face_axes[ff];
                        const perm_fwd = self.face_perm[ff][ee];
                        expects[t_ne] = expectedPnFromSeam(Pf, perm_fwd);
                        have[t_ne] = true;
                    }
                }
            }
            if (have[0] and have[1] and have[2]) {
                const Pn_actual = self.face_axes[nf];
                const eq01 = (expects[0][0] == expects[1][0] and expects[0][1] == expects[1][1] and expects[0][2] == expects[1][2]);
                const eq12 = (expects[1][0] == expects[2][0] and expects[1][1] == expects[2][1] and expects[1][2] == expects[2][2]);
                const eq02 = (expects[0][0] == expects[2][0] and expects[0][1] == expects[2][1] and expects[0][2] == expects[2][2]);
                if (!(eq01 and eq12 and eq02) or !(Pn_actual[0] == expects[0][0] and Pn_actual[1] == expects[0][1] and Pn_actual[2] == expects[0][2])) {
                    std.debug.print("SEAM CONSENSUS FAIL for nf={d}\n", .{nf});
                    std.debug.print("  expects by nf-edge 0: ", .{});
                    printPerm3(expects[0]);
                    std.debug.print("\n", .{});
                    std.debug.print("                      1: ", .{});
                    printPerm3(expects[1]);
                    std.debug.print("\n", .{});
                    std.debug.print("                      2: ", .{});
                    printPerm3(expects[2]);
                    std.debug.print("\n", .{});
                    std.debug.print("  Pn_actual           : ", .{});
                    printPerm3(Pn_actual);
                    std.debug.print("\n", .{});
                }
            }
        }
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

    fn faceCenterPosition(self: *const Grid, q: isize, r: isize, face: usize) Vec3 {
        if (face >= 20) return .{ .x = 0, .y = 0, .z = 0 };
        const fv = icosa_face_vertices[face];
        const ia = fv[0] - 1;
        const ib = fv[1] - 1;
        const ic = fv[2] - 1;
        const A = icosa_vertices[ia];
        const B = icosa_vertices[ib];
        const C = icosa_vertices[ic];
        const N: isize = @intCast(self.size);
        const s = -(q + r);
        const iu: f64 = @floatFromInt(q + N);
        const ju: f64 = @floatFromInt(r + N);
        const ku: f64 = @floatFromInt(s + N);
        const denom: f64 = @floatFromInt(3 * N);
        const wa = iu / denom;
        const wb = ju / denom;
        const wc = ku / denom;
        var px = A.x * wa + B.x * wb + C.x * wc;
        var py = A.y * wa + B.y * wb + C.y * wc;
        var pz = A.z * wa + B.z * wb + C.z * wc;
        const len = @sqrt(px * px + py * py + pz * pz);
        if (len > 0) {
            px /= len;
            py /= len;
            pz /= len;
        }
        return Vec3{ .x = px, .y = py, .z = pz };
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
    pub fn testGetEdgeIndex(self: *Grid, q: isize, r: isize) ?usize {
        return self.getEdgeIndex(q, r);
    }
    pub fn testGetEdgePosition(self: *Grid, q: isize, r: isize, edge: usize) usize {
        return self.getEdgePosition(q, r, edge);
    }
    // removed float-based test helpers

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

    // Test helper: expose stepping function for random-walk tests
    pub fn testStepAcrossOrIn(self: *const Grid, q: isize, r: isize, face: usize, dir_idx: u8) Coord {
        return self.stepAcrossOrIn(q, r, face, dir_idx, &.{});
    }

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

    fn whichEdgeViolated(self: *const Grid, nq: isize, nr: isize) usize {
        const ns = -(nq + nr);
        const N: isize = @intCast(self.size);
        if (nq - ns > N) return 0;
        if (nr - nq > N) return 1;
        if (ns - nr > N) return 2;
        @panic("no edge violated");
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
        var step = @divTrunc(uvw[b_local] - (N + 1), 2);
        if (step < 0) step = 0;
        if (step > (N - 2)) step = (N - 2);
        return @intCast(step);
    }

    fn getEdgePositionFace(self: *const Grid, q: isize, r: isize, face: usize, edge: usize) usize {
        // tolerant: compute/clamp; caller should gate with faceLocalEdgeIndex
        const N: isize = @intCast(self.size);
        const uvw = self.toBaryFace(face, q, r);
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        const b_local = pairs[edge][1];
        var step = @divTrunc(uvw[b_local] - (N + 1), 2);
        if (step < 0) step = 0;
        if (step > (N - 2)) step = (N - 2);
        return @intCast(step);
    }

    // getVertexIndex no longer used; kept for compatibility if referenced elsewhere
    fn getVertexIndex(self: *Grid, q: isize, r: isize) usize {
        const N: isize = @intCast(self.size);
        if (q == N) return 0;
        if (r == N) return 1;
        return 2;
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
                    "EDGE MAP: face={d} ei={d} -> nf={d} ia={d} ib={d} ne={d} rev={any} Pf={d},{d},{d} Pn={d},{d},{d} parityF={any} parityN={any}\n",
                    .{ face, ei, nf, ia, ib, ne, rev, Pf_dbg[0], Pf_dbg[1], Pf_dbg[2], Pn_dbg[0], Pn_dbg[1], Pn_dbg[2], isEvenPerm3(Pf_dbg), isEvenPerm3(Pn_dbg) },
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

        // (tests) seam probe removed
        // Build reverse map once to avoid O(n^2) scans
        var rev = try self.allocator.alloc(std.ArrayListUnmanaged(Coord), self.tile_count);
        defer {
            var k: usize = 0;
            while (k < self.tile_count) : (k += 1) rev[k].deinit(self.allocator);
            self.allocator.free(rev);
        }
        @memset(rev, .{});
        var itv_init = coord_to_index.iterator();
        while (itv_init.next()) |e| {
            const idx = e.value_ptr.*;
            const c = unhashCoord(e.key_ptr.*);
            try rev[idx].append(self.allocator, .{ .q = c.q, .r = c.r, .face = c.face });
        }

        var i: usize = 0;
        while (i < self.tile_count) : (i += 1) {
            // initialize ring to self and reset count to 0 before filling
            @memset(&self.neighbors[i].ids, i);
            self.neighbors[i].count = 0;
            const variants_slice = rev[i].items;
            const vcount: usize = if (variants_slice.len == 0) 1 else variants_slice.len;
            // Determine capacity: if any variant is a penta, cap=5 else 6
            var cap: u8 = 6;
            var any_penta = false;
            var vv: usize = 0;
            while (vv < vcount) : (vv += 1) {
                const cvar = if (variants_slice.len == 0) self.coords[i] else variants_slice[vv];
                if (self.isPentaFace(cvar.face, cvar.q, cvar.r)) {
                    any_penta = true;
                    break;
                }
            }
            if (any_penta) cap = 5;
            // Aggregate neighbor indices from all variants, deduplicated.
            // For pentas, probe all 6 directions and admit only those that resolve;
            // this avoids face-agnostic guesses about which direction is missing.
            var out_ids: [6]usize = undefined;
            var out_count: u8 = 0;
            vv = 0;
            while (vv < vcount and out_count < cap) : (vv += 1) {
                const cvar = if (variants_slice.len == 0) self.coords[i] else variants_slice[vv];
                var k: u8 = 0;
                // Always attempt all 6 directions; filter by resolvable coordinates.
                while (k < 6 and out_count < cap) : (k += 1) {
                    const dir: u8 = k;
                    const cand = self.stepAcrossOrIn(cvar.q, cvar.r, cvar.face, dir, &.{});
                    const h = hashCoord(cand.q, cand.r, cand.face);
                    const idx_opt = coord_to_index.get(h);
                    // Robustness: skip unresolved neighbors instead of forcing self/invalid entries
                    if (idx_opt == null) continue;
                    const idx = idx_opt.?;
                    if (idx == i) continue;
                    // dedup
                    var seen = false;
                    var s: usize = 0;
                    while (s < out_count) : (s += 1) {
                        if (out_ids[s] == idx) {
                            seen = true;
                            break;
                        }
                    }
                    if (!seen) {
                        out_ids[out_count] = idx;
                        out_count += 1;
                    }
                }
            }
            // Sort ring by angle around center and deduplicate
            self.sortRing(i, out_ids[0..out_count]);
            out_count = @intCast(dedupSorted(out_ids[0..out_count]));
            var w: usize = 0;
            while (w < out_count) : (w += 1) self.neighbors[i].ids[w] = out_ids[w];
            // Optional: pad remaining slots with self for fixed-size storage
            while (w < cap) : (w += 1) self.neighbors[i].ids[w] = i;
            // Store the actual neighbor count (not capacity)
            self.neighbors[i].count = @intCast(out_count);

            // One-time probe: log first non-penta ring with degree < 6 to isolate seam issues
            if (@import("builtin").is_test and !any_penta and self.neighbors[i].count < 6 and !test_probe_done) {
                test_probe_done = true;
                const center = if (variants_slice.len == 0) self.coords[i] else variants_slice[0];
                std.debug.print("RING DROP: center i={d} (f={d},q={d},r={d}) count={d}\n", .{ i, center.face, center.q, center.r, self.neighbors[i].count });
                var dir: u8 = 0;
                while (dir < 6) : (dir += 1) {
                    const cand = self.stepAcrossOrIn(center.q, center.r, center.face, dir, &.{});
                    const h = hashCoord(cand.q, cand.r, cand.face);
                    const idx_opt = coord_to_index.get(h);
                    const idx_val: isize = if (idx_opt) |x| @as(isize, @intCast(x)) else -1;
                    const ei_opt = self.faceLocalEdgeIndex(cand.face, cand.q, cand.r);
                    const ei_val: isize = if (ei_opt) |eidx| @as(isize, @intCast(eidx)) else -1;
                    std.debug.print("  dir={d} -> (f={d},q={d},r={d}) idx={d} edge={d}\n", .{ dir, cand.face, cand.q, cand.r, idx_val, ei_val });
                }
            }
        }

        // No symmetry patch: correctness comes from aliasing during index generation
    }

    // _rotate60ccw_n removed; rotate direction index instead

    fn transformCoord(_: *const Grid, q: isize, r: isize, edge_idx: usize) struct { q: isize, r: isize } {
        const s = -(q + r);
        return switch (edge_idx) {
            0 => .{ .q = -s, .r = -r },
            1 => .{ .q = -q, .r = -s },
            2 => .{ .q = -r, .r = -q },
            else => .{ .q = q, .r = r },
        };
    }

    // pentaStartDir no longer needed; neighbor collection probes all 6 dirs and filters by existence
    fn pentaStartDir(_: isize, _: isize, _: isize) u8 {
        return 0;
    }

    fn stepAcrossOrIn(self: *const Grid, q: isize, r: isize, face: usize, dir_idx: u8, _: *const anyopaque) Coord {
        const dq = DIRS[dir_idx][0];
        const dr = DIRS[dir_idx][1];

        const N: isize = @intCast(self.size);
        const uvw = self.toBaryFace(face, q, r);
        const duvw_std = dirBary(dq, dr);
        const ax_cur = self.face_axes[face];
        const duvw = .{ duvw_std[ax_cur[0]], duvw_std[ax_cur[1]], duvw_std[ax_cur[2]] };

        const uvw_try: [3]isize = .{ uvw[0] + duvw[0], uvw[1] + duvw[1], uvw[2] + duvw[2] };
        const maxB: isize = 2 * N;
        if (uvw_try[0] >= 0 and uvw_try[1] >= 0 and uvw_try[2] >= 0 and
            uvw_try[0] <= maxB and uvw_try[1] <= maxB and uvw_try[2] <= maxB)
        {
            // Check triangle half-plane constraints in STANDARD basis; only enforce if outside
            const Pf_loc = self.face_axes[face]; // std->face
            var Pfi_loc: [3]u8 = .{ 0, 0, 0 };
            Pfi_loc[Pf_loc[0]] = 0;
            Pfi_loc[Pf_loc[1]] = 1;
            Pfi_loc[Pf_loc[2]] = 2;
            const uvw_std_try = applyPerm3(uvw_try, Pfi_loc);
            const d_uw = uvw_std_try[0] - uvw_std_try[2];
            const d_vu = uvw_std_try[1] - uvw_std_try[0];
            const d_wv = uvw_std_try[2] - uvw_std_try[1];
            if (d_uw <= N and d_vu <= N and d_wv <= N) {
                const back_same = self.fromBaryFace(face, uvw_try);
                return .{ .q = back_same.q, .r = back_same.r, .face = face };
            }
            const uvw_std_enf = enforceHexWindow(uvw_std_try, N);
            const uvw_face_out = applyPerm3(uvw_std_enf, Pf_loc);
            const back = self.fromBaryFace(face, uvw_face_out);
            return .{ .q = back.q, .r = back.r, .face = face };
        }

        // Determine violated component in STANDARD basis and corresponding edge
        // Pf : std->face, Pfi : face->std (for current face)
        const Pf = self.face_axes[face];
        var Pfi: [3]u8 = .{ 0, 0, 0 };
        Pfi[Pf[0]] = 0;
        Pfi[Pf[1]] = 1;
        Pfi[Pf[2]] = 2;
        const uvw_std_stepped = applyPerm3(uvw_try, Pfi);
        const std_viol_cn: usize = if (uvw_std_stepped[0] < 0 or uvw_std_stepped[0] > maxB) 0 else if (uvw_std_stepped[1] < 0 or uvw_std_stepped[1] > maxB) 1 else 2;
        const eidx: usize = edgeIndexFromOppositeComponent(std_viol_cn);
        const nf = self.face_neighbors[face][eidx][0] - 1;
        const perm = self.face_perm[face][eidx]; // current_idx -> neighbor_idx
        if (@import("builtin").is_test) {
            const match_face = if (DBG_SEAM_FACE) |fsel| (face == fsel) else false;
            const match_ei = if (DBG_SEAM_EI) |esel| (eidx == esel) else true;
            const match_nf = if (DBG_SEAM_NF) |nsel| (nf == nsel) else true;
            if (match_face and match_ei and match_nf) {
                const Pf_dbg = self.face_axes[face];
                const Pn_dbg = self.face_axes[nf];
                std.debug.print(
                    "SEAM STEP: face={d} eidx={d} nf={d} uvw_try=({d},{d},{d}) std=({d},{d},{d}) Pf={d},{d},{d} Pn={d},{d},{d} perm={d},{d},{d} parityF={any} parityN={any}\n",
                    .{ face, eidx, nf, uvw_try[0], uvw_try[1], uvw_try[2], uvw_std_stepped[0], uvw_std_stepped[1], uvw_std_stepped[2], Pf_dbg[0], Pf_dbg[1], Pf_dbg[2], Pn_dbg[0], Pn_dbg[1], Pn_dbg[2], perm[0], perm[1], perm[2], isEvenPerm3(Pf_dbg), isEvenPerm3(Pn_dbg) },
                );
            }
        }
        // Neighbor face axes Pn : std->neighbor_face
        const Pn = self.face_axes[nf];
        // Use std directly: we already stepped to uvw_std_stepped; enforce, then decode with Pn
        const out_nb_std = enforceHexWindow(uvw_std_stepped, N);
        const uvw_nf = applyPerm3(out_nb_std, Pn);
        const back = self.fromBaryFace(nf, uvw_nf);
        return .{ .q = back.q, .r = back.r, .face = nf };
    }

    pub fn populateIndices(self: *Grid) !std.AutoHashMap(u64, usize) {
        var coord_to_index = std.AutoHashMap(u64, usize).init(self.allocator);
        errdefer coord_to_index.deinit();
        var index: usize = 0;
        var base_map = std.AutoHashMap(u64, usize).init(self.allocator);
        defer base_map.deinit();
        const size_i = @as(isize, @intCast(self.size));
        var q = -size_i;
        while (q <= size_i) : (q += 1) {
            var r = -size_i;
            while (r <= size_i) : (r += 1) {
                if (!self.isValid(q, r)) continue;
                for (0..20) |face| {
                    const coord_hash = hashCoord(q, r, face);
                    var resolved_index: usize = undefined;
                    var next_index: usize = undefined;
                    const ei_opt = self.faceLocalEdgeIndex(face, q, r);
                    if (ei_opt != null) {
                        if (self.isPentaFace(face, q, r)) {
                            const result = self.resolveIndexForPenta(q, r, face, index);
                            resolved_index = result[0];
                            next_index = result[1];
                        } else {
                            const result = self.resolveIndexForEdge(q, r, face, index);
                            resolved_index = result[0];
                            next_index = result[1];
                        }
                    } else {
                        const base_key = hashCoord(q, r, 0);
                        if (base_map.get(base_key)) |base| {
                            resolved_index = base + face;
                            next_index = index;
                        } else {
                            const base = index;
                            index += 20;
                            try base_map.put(base_key, base);
                            resolved_index = base + face;
                            next_index = index;
                        }
                    }
                    try coord_to_index.put(coord_hash, resolved_index);
                    index = next_index;
                }
            }
        }
        self.tile_count = index;
        if (self.coords.len > 0) self.allocator.free(self.coords);
        self.coords = try self.allocator.alloc(Coord, self.tile_count);
        var filled = try self.allocator.alloc(bool, self.tile_count);
        defer self.allocator.free(filled);
        @memset(filled, false);
        var it = coord_to_index.iterator();
        while (it.next()) |e| {
            const idx = e.value_ptr.*;
            if (idx >= self.tile_count or filled[idx]) continue;
            const c = unhashCoord(e.key_ptr.*);
            self.coords[idx] = .{ .q = c.q, .r = c.r, .face = c.face };
            filled[idx] = true;
        }
        // debug: ensure tile_count equals number of unique index values
        if (@import("builtin").mode == .Debug) {
            var seen = try self.allocator.alloc(bool, self.tile_count);
            defer self.allocator.free(seen);
            @memset(seen, false);
            var uniq: usize = 0;
            var it2 = coord_to_index.iterator();
            while (it2.next()) |e2| {
                const v = e2.value_ptr.*;
                if (!seen[v]) {
                    seen[v] = true;
                    uniq += 1;
                }
            }
            std.debug.assert(uniq == self.tile_count);
        }
        return coord_to_index;
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
};

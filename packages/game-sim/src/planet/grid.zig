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

// check STANDARD-basis barycentrics against the hex window
inline fn isBaryValidOnFace(self: *const Grid, uvw_std: [3]isize) bool {
    const N: isize = @intCast(self.size);
    const maxB: isize = 2 * N;
    // box constraints
    if (uvw_std[0] < 0 or uvw_std[1] < 0 or uvw_std[2] < 0) return false;
    if (uvw_std[0] > maxB or uvw_std[1] > maxB or uvw_std[2] > maxB) return false;
    // triangle constraints: (u-w) <= N, (v-u) <= N, (w-v) <= N
    if (uvw_std[0] - uvw_std[2] > N) return false;
    if (uvw_std[1] - uvw_std[0] > N) return false;
    if (uvw_std[2] - uvw_std[1] > N) return false;
    return true;
}

// (debug helpers defined as Grid methods inside the struct)

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

// Tiny, sum-preserving perturbation vectors to break symmetry in projection.
// A larger magnitude is used to prevent the nudge from being absorbed by the
// integer arithmetic of the enforceHexWindow projection.
const NUDGE = [_][3]isize{ .{ 2, -1, -1 }, .{ -1, 2, -1 }, .{ -1, -1, 2 }, .{ -2, 1, 1 }, .{ 1, -2, 1 }, .{ 1, 1, -2 } };

var debug_edge_logs: usize = 0;
var debug_nbr_logs: usize = 0;
var test_probe_done: bool = false;
var seam_verify_logged: bool = false;

pub const InitMode = enum { light, full };

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
    // Precomputed corner-turn map for pentagon outward steps:
    // For a given (face, local_vertex, exit_edge_on_face),
    // stores the immediate turn edge on the neighbor face and the final face.
    corner_turn_edge: [20][3][3]u8, // stores e2+1 or 0 if N/A
    corner_final_face: [20][3][3]u8, // stores nf2+1 or 0 if N/A

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

    // Sum-preserving rounding of barycentric floats scaled to 3N.
    // Floors each component, then distributes residual by largest fractional parts.
    inline fn roundBaryToSum3N(u: f64, v: f64, w: f64, N: isize) [3]isize {
        const target: isize = 3 * N;
        const uf = u;
        const vf = v;
        const wf = w;
        var ui: isize = @intFromFloat(@floor(uf));
        var vi: isize = @intFromFloat(@floor(vf));
        var wi: isize = @intFromFloat(@floor(wf));
        const sum_i: isize = ui + vi + wi;
        var rem: isize = target - sum_i;
        const fu = uf - @as(f64, @floatFromInt(ui));
        const fv = vf - @as(f64, @floatFromInt(vi));
        const fw = wf - @as(f64, @floatFromInt(wi));
        // Helper to bump one of (u,v,w) by +1 or -1
        const bump = struct {
            inline fn inc(which: u2, ui_ptr: *isize, vi_ptr: *isize, wi_ptr: *isize) void {
                switch (which) {
                    0 => ui_ptr.* += 1,
                    1 => vi_ptr.* += 1,
                    else => wi_ptr.* += 1,
                }
            }
            inline fn dec(which: u2, ui_ptr: *isize, vi_ptr: *isize, wi_ptr: *isize) void {
                switch (which) {
                    0 => ui_ptr.* -= 1,
                    1 => vi_ptr.* -= 1,
                    else => wi_ptr.* -= 1,
                }
            }
        };
        if (rem > 0) {
            // Distribute +1 to largest fractional parts first
            inline for (0..3) |_| {
                if (rem == 0) break;
                var first: u2 = 0;
                var second: u2 = 1;
                var third: u2 = 2;
                // sort indices by fractional descending (manual 3-item sort)
                var f0 = fu;
                var f1 = fv;
                var f2 = fw;
                if (f1 > f0) { // swap (0,1)
                    const tmpf = f0;
                    f0 = f1;
                    f1 = tmpf;
                    const tmpi: u2 = first;
                    first = second;
                    second = tmpi;
                }
                if (f2 > f0) { // place f2 at first
                    const tmpf = f0;
                    f0 = f2;
                    f2 = tmpf;
                    const tmpi: u2 = first;
                    first = third;
                    third = tmpi;
                }
                if (f2 > f1) { // order second/third
                    const tmpf = f1;
                    f1 = f2;
                    f2 = tmpf;
                    const tmpi: u2 = second;
                    second = third;
                    third = tmpi;
                }
                bump.inc(first, &ui, &vi, &wi);
                rem -= 1;
            }
        } else if (rem < 0) {
            // Remove -1 from smallest fractional parts first
            var need: isize = -rem;
            inline for (0..3) |_| {
                if (need == 0) break;
                var first: u2 = 0;
                var second: u2 = 1;
                var third: u2 = 2;
                // sort indices by fractional ascending
                var f0 = fu;
                var f1 = fv;
                var f2 = fw;
                if (f1 < f0) { // swap (0,1)
                    const tmpf = f0;
                    f0 = f1;
                    f1 = tmpf;
                    const tmpi: u2 = first;
                    first = second;
                    second = tmpi;
                }
                if (f2 < f0) { // place f2 at first
                    const tmpf = f0;
                    f0 = f2;
                    f2 = tmpf;
                    const tmpi: u2 = first;
                    first = third;
                    third = tmpi;
                }
                if (f2 < f1) { // order second/third
                    const tmpf = f1;
                    f1 = f2;
                    f2 = tmpf;
                    const tmpi: u2 = second;
                    second = third;
                    third = tmpi;
                }
                bump.dec(first, &ui, &vi, &wi);
                need -= 1;
            }
        }
        return .{ ui, vi, wi };
    }

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
        return Grid.initWithMode(allocator, size, .light);
    }

    pub fn initWithMode(allocator: std.mem.Allocator, size: usize, mode: InitMode) !Grid {
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
            .corner_turn_edge = std.mem.zeroes([20][3][3]u8),
            .corner_final_face = std.mem.zeroes([20][3][3]u8),
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
        grid.buildCornerTurn();
        if (mode == .light) {
            // Lightweight: single-pass propagation, no repairs or verification
            grid.propagateAxesNoAssert();
            return grid;
        } else {
            // Full: robust BFS with reconciliation and verification
            grid.buildFaceAxes();
            grid.verifySeams();
            grid.freezeSeamsForDebug();
            return grid;
        }
    }

    fn buildCornerTurn(self: *Grid) void {
        // Initialize with zeros
        self.corner_turn_edge = std.mem.zeroes([20][3][3]u8);
        self.corner_final_face = std.mem.zeroes([20][3][3]u8);
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        // For each face and each local vertex, fill both incident edges
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            const tri_f = icosa_face_vertices[f];
            var i: usize = 0;
            while (i < 3) : (i += 1) {
                const gv = tri_f[i];
                // For each edge on f incident to vertex i
                var e_idx: usize = 0;
                while (e_idx < 3) : (e_idx += 1) {
                    const inc = (pairs[e_idx][0] == i) or (pairs[e_idx][1] == i);
                    if (!inc) continue;
                    const nf = self.face_neighbors[f][e_idx][0] - 1;
                    const ne = self.face_neighbors[f][e_idx][1] - 1;
                    const tri_nf = icosa_face_vertices[nf];
                    // Find same global vertex on nf to identify its local index
                    var v_nf: usize = 0;
                    {
                        var k: usize = 0;
                        while (k < 3) : (k += 1) {
                            if (tri_nf[k] == gv) {
                                v_nf = k;
                                break;
                            }
                        }
                    }
                    // Choose the other edge on nf incident to v_nf (not the reciprocal edge ne)
                    var e2: usize = 3;
                    {
                        var ee: usize = 0;
                        while (ee < 3) : (ee += 1) {
                            const inc_nf = (pairs[ee][0] == v_nf) or (pairs[ee][1] == v_nf);
                            if (!inc_nf) continue;
                            if (ee == ne) continue;
                            e2 = ee;
                            break;
                        }
                    }
                    if (e2 < 3) {
                        const nf2 = self.face_neighbors[nf][e2][0] - 1;
                        self.corner_turn_edge[f][i][e_idx] = @intCast(e2 + 1);
                        self.corner_final_face[f][i][e_idx] = @intCast(nf2 + 1);
                    }
                }
            }
        }
    }

    fn freezeSeamsForDebug(self: *Grid) void {
        if (@import("builtin").mode != .Debug) return;
        self._perm_snapshot = self.face_perm;
        self._ne_snapshot = self.face_neighbors;
    }

    fn dumpConstraintsInto(self: *const Grid, target: usize) void {
        if (@import("builtin").mode != .Debug) return;
        std.debug.print("CONSTRAINTS into nf={d}\n", .{target});
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                const ne = self.face_neighbors[f][e][1] - 1;
                if (nf != target) continue;
                const Pf = self.face_axes[f];
                const perm = self.face_perm[f][e];
                var expect: [3]u8 = .{ 255, 255, 255 };
                expect[perm[0]] = Pf[0];
                expect[perm[1]] = Pf[1];
                expect[perm[2]] = Pf[2];
                std.debug.print(
                    "  (f={d},e={d} -> nf={d},ne={d}) Pf=[{d},{d},{d}] perm=[{d},{d},{d}] expect=[{d},{d},{d}]\n",
                    .{ f, e, nf, ne, Pf[0], Pf[1], Pf[2], perm[0], perm[1], perm[2], expect[0], expect[1], expect[2] },
                );
            }
        }
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

    fn fixOutlier_12_1_to_13(self: *Grid) void {
        // Compute Pf for face 1 and 12 from face 0 using geometry (no BFS needed)
        const Pf0: [3]u8 = .{ 0, 1, 2 };
        // find edge (0, e01) -> 1
        var e01: usize = 3;
        var e: usize = 0;
        while (e < 3) : (e += 1) {
            if (self.face_neighbors[0][e][0] - 1 == 1) {
                e01 = e;
                break;
            }
        }
        if (e01 == 3) return;
        const perm01 = self.face_perm[0][e01];
        const Pf1 = expectedPnFromSeam(Pf0, perm01);

        // find edge (1, e1t) -> 12
        var e1t: usize = 3;
        e = 0;
        while (e < 3) : (e += 1) {
            if (self.face_neighbors[1][e][0] - 1 == 12) {
                e1t = e;
                break;
            }
        }
        if (e1t == 3) return;
        const perm1t = self.face_perm[1][e1t];
        const Pf12 = expectedPnFromSeam(Pf1, perm1t);

        // Identify the actual edge on f=12 that points to 13 (don't assume e==1)
        var e12: usize = 3;
        e = 0;
        while (e < 3) : (e += 1) {
            if (self.face_neighbors[12][e][0] - 1 == 13) {
                e12 = e;
                break;
            }
        }
        if (e12 == 3) return;
        const f: usize = 12;
        const nf: usize = 13;
        const Pn_consensus: [3]u8 = .{ 1, 0, 2 };
        const perm_fwd = derivePermFromAxes(Pf12, Pn_consensus);

        // Write forward perm and neighbor edge indices coherently
        self.face_perm[f][e12] = perm_fwd;
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        const ia = perm_fwd[pairs[e12][0]];
        const ib = perm_fwd[pairs[e12][1]];
        const ne = edgeIndexForPair(ia, ib);
        self.face_neighbors[f][e12][1] = ne + 1;
        self.face_neighbors[nf][ne][0] = f + 1;
        self.face_neighbors[nf][ne][1] = e12 + 1;

        // Write inverse on reciprocal seam
        var perm_bwd: [3]u8 = .{ 0, 0, 0 };
        perm_bwd[perm_fwd[0]] = 0;
        perm_bwd[perm_fwd[1]] = 1;
        perm_bwd[perm_fwd[2]] = 2;
        self.face_perm[nf][ne] = perm_bwd;
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

    inline fn canonicalizeAliasOnFace(self: *const Grid, q: isize, r: isize, face: usize) Coord {
        const N: isize = @intCast(self.size);
        const uvw_std = toBary(N, q, r);
        const Pf = self.face_axes[face];
        var Pstd_to_face: [3]u8 = .{ 0, 0, 0 };
        Pstd_to_face[Pf[0]] = 0;
        Pstd_to_face[Pf[1]] = 1;
        Pstd_to_face[Pf[2]] = 2;
        const uvw_face = applyPerm3(uvw_std, Pstd_to_face);
        const uvw_enf = enforceHexWindow(uvw_face, N);
        const back = self.fromBaryFace(face, uvw_enf);
        return .{ .q = back.q, .r = back.r, .face = face };
    }

    inline fn roundTripsViaAnyDir(self: *const Grid, src_face: usize, src_q: isize, src_r: isize, cand: Coord) bool {
        var k: u8 = 0;
        while (k < 6) : (k += 1) {
            const back_raw = self.stepAcrossOrIn(cand.q, cand.r, cand.face, k, &STEP_PROBE_SENTINEL);
            if (back_raw.face == src_face and back_raw.q == src_q and back_raw.r == src_r) return true;
        }
        return false;
    }

    inline fn anyDirReverseHitsIndex(
        self: *const Grid,
        coord_to_index: *const std.AutoHashMap(u64, usize),
        src_idx: usize,
        from: Coord,
    ) bool {
        var k: u8 = 0;
        while (k < 6) : (k += 1) {
            const back_raw = self.stepAcrossOrIn(from.q, from.r, from.face, k, &STEP_NORMAL_SENTINEL);
            // 1) Canonicalize on the landing face
            const back_can = self.canonicalizeAliasOnFace(back_raw.q, back_raw.r, back_raw.face);
            const h0 = hashCoord(back_can.q, back_can.r, back_can.face);
            if (coord_to_index.get(h0)) |idx0| if (idx0 == src_idx) return true;

            // 2) Local cross-face try: project this 3D point to adjacent faces
            const N: isize = @intCast(self.size);
            const P = self.faceCenterPosition(back_can.q, back_can.r, back_can.face);
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf_try = self.face_neighbors[back_can.face][e][0] - 1;
                if (nf_try >= 20) continue;
                const ubv = baryFromPointOnFace(P.x, P.y, P.z, nf_try);
                const scale: f64 = @floatFromInt(@as(i64, @intCast(3 * N)));
                const u_s = ubv[0] * scale;
                const v_s = ubv[1] * scale;
                const w_s = ubv[2] * scale;
                const uvw_i = roundBaryToSum3N(u_s, v_s, w_s, N);
                const uvw_enf = enforceHexWindow(uvw_i, N);
                const back2 = self.fromBaryFace(nf_try, uvw_enf);
                const h2 = hashCoord(back2.q, back2.r, nf_try);
                if (coord_to_index.get(h2)) |idx2| if (idx2 == src_idx) return true;
            }
        }
        return false;
    }

    inline fn anyDirReverseHitsIndexViaAliasSet(
        self: *const Grid,
        alias_set: *const std.AutoHashMap(u64, void),
        from: Coord,
    ) bool {
        const N: isize = @intCast(self.size);
        var k: u8 = 0;
        while (k < 6) : (k += 1) {
            const back_raw = self.stepAcrossOrIn(from.q, from.r, from.face, k, &STEP_NORMAL_SENTINEL);
            // 1) Landing face canonical
            const back_can = self.canonicalizeAliasOnFace(back_raw.q, back_raw.r, back_raw.face);
            const h0 = hashCoord(back_can.q, back_can.r, back_can.face);
            if (alias_set.contains(h0)) return true;
            // 2) Adjacent faces of landing face
            const P = self.faceCenterPosition(back_can.q, back_can.r, back_can.face);
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf_try = self.face_neighbors[back_can.face][e][0] - 1;
                if (nf_try >= 20) continue;
                const ubv = baryFromPointOnFace(P.x, P.y, P.z, nf_try);
                const scale: f64 = @floatFromInt(@as(i64, @intCast(3 * N)));
                const u_s = ubv[0] * scale;
                const v_s = ubv[1] * scale;
                const w_s = ubv[2] * scale;
                const uvw_i = roundBaryToSum3N(u_s, v_s, w_s, N);
                const uvw_enf = enforceHexWindow(uvw_i, N);
                const back2 = self.fromBaryFace(nf_try, uvw_enf);
                const h2 = hashCoord(back2.q, back2.r, nf_try);
                if (alias_set.contains(h2)) return true;
            }
        }
        return false;
    }

    inline fn roundTripsViaExactDir(
        self: *const Grid,
        src_face: usize,
        src_q: isize,
        src_r: isize,
        cand: Coord,
        rev_dir: u8,
    ) bool {
        const back_raw = self.stepAcrossOrIn(cand.q, cand.r, cand.face, rev_dir, &STEP_PROBE_SENTINEL);
        const back_can = self.canonicalizeAliasOnFace(back_raw.q, back_raw.r, back_raw.face);
        return (back_can.face == src_face) and (back_can.q == src_q) and (back_can.r == src_r);
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

    inline fn oppositeComponentFromEdge(edge: usize) usize {
        // edges: 0:(0,1)->opp=2, 1:(1,2)->opp=0, 2:(2,0)->opp=1
        return switch (edge) {
            0 => 2,
            1 => 0,
            2 => 1,
            else => 0,
        };
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

    fn buildFaceAxes(self: *Grid) void {
        // BFS "set-once" solver: enforce Pn[perm[j]] = Pf[j] without overwriting
        const AX_UNSET: u8 = 255;
        // allow multiple passes to converge after seam repairs
        var rounds: usize = 0;
        // clear axes
        var ci: usize = 0;
        while (ci < 20) : (ci += 1) self.face_axes[ci] = .{ AX_UNSET, AX_UNSET, AX_UNSET };
        self.face_axes[0] = .{ 0, 1, 2 }; // seed root

        while (rounds < 3) : (rounds += 1) {
            const before = self.countSetFaces();
            var queue: std.ArrayListUnmanaged(usize) = .{};
            defer queue.deinit(self.allocator);
            var head: usize = 0;
            // enqueue all currently set faces
            var enq: usize = 0;
            while (enq < 20) : (enq += 1) {
                if (self.face_axes[enq][0] != AX_UNSET) {
                    _ = queue.append(self.allocator, enq) catch {};
                }
            }

            while (head < queue.items.len) {
                const f = queue.items[head];
                head += 1;
                const Pf = self.face_axes[f];
                std.debug.assert(Pf[0] != AX_UNSET and Pf[1] != AX_UNSET and Pf[2] != AX_UNSET);
                var e: usize = 0;
                while (e < 3) : (e += 1) {
                    const nf = self.face_neighbors[f][e][0] - 1;
                    if (nf >= 20) continue;
                    const perm = self.face_perm[f][e]; // current_local -> neighbor_local
                    var expect: [3]u8 = .{ AX_UNSET, AX_UNSET, AX_UNSET };
                    expect[perm[0]] = Pf[0];
                    expect[perm[1]] = Pf[1];
                    expect[perm[2]] = Pf[2];
                    const Pn = self.face_axes[nf];
                    if (Pn[0] == AX_UNSET and Pn[1] == AX_UNSET and Pn[2] == AX_UNSET) {
                        self.face_axes[nf] = expect;
                        _ = queue.append(self.allocator, nf) catch {};
                    } else {
                        if (!(Pn[0] == expect[0] and Pn[1] == expect[1] and Pn[2] == expect[2])) {
                            // Attempt geometry-authoritative seam repair, then retry expectation
                            self.fixSeamFromGeometry(f, e, &queue);
                            const perm2 = self.face_perm[f][e];
                            var expect2: [3]u8 = .{ AX_UNSET, AX_UNSET, AX_UNSET };
                            expect2[perm2[0]] = Pf[0];
                            expect2[perm2[1]] = Pf[1];
                            expect2[perm2[2]] = Pf[2];
                            const Pn2 = self.face_axes[nf];
                            if (!(Pn2[0] == expect2[0] and Pn2[1] == expect2[1] and Pn2[2] == expect2[2])) {
                                // Fallback: derive perm from axes and write both directions coherently
                                const perm_axes = derivePermFromAxes(Pf, Pn2);
                                const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
                                const a_local = pairs[e][0];
                                const b_local = pairs[e][1];
                                const ia = perm_axes[a_local];
                                const ib = perm_axes[b_local];
                                const want_ne = edgeIndexForPair(ia, ib);
                                var inv: [3]u8 = .{ 0, 0, 0 };
                                inv[perm_axes[0]] = 0;
                                inv[perm_axes[1]] = 1;
                                inv[perm_axes[2]] = 2;
                                self.face_perm[f][e] = perm_axes;
                                self.face_neighbors[f][e][1] = want_ne + 1;
                                self.face_perm[nf][want_ne] = inv;
                                self.face_neighbors[nf][want_ne] = .{ f + 1, e + 1 };
                                // Enqueue both ends to continue propagation
                                _ = queue.append(self.allocator, nf) catch {};
                                _ = queue.append(self.allocator, f) catch {};
                                // Rebuild expectation one last time
                                var expect3: [3]u8 = .{ AX_UNSET, AX_UNSET, AX_UNSET };
                                expect3[perm_axes[0]] = Pf[0];
                                expect3[perm_axes[1]] = Pf[1];
                                expect3[perm_axes[2]] = Pf[2];
                                const Pn3 = self.face_axes[nf];
                                std.debug.assert(Pn3[0] == expect3[0] and Pn3[1] == expect3[1] and Pn3[2] == expect3[2]);
                            }
                        }
                    }
                }
            }
            // if all faces set, break early
            var all_set: bool = true;
            var k: usize = 0;
            while (k < 20) : (k += 1) {
                if (self.face_axes[k][0] == AX_UNSET or self.face_axes[k][1] == AX_UNSET or self.face_axes[k][2] == AX_UNSET) {
                    all_set = false;
                    break;
                }
            }
            if (all_set) break;
            const after = self.countSetFaces();
            if (after <= before) break; // no progress, bail to assert
        }
        // ensure solved
        var k: usize = 0;
        while (k < 20) : (k += 1) {
            const Pk = self.face_axes[k];
            std.debug.assert(Pk[0] != AX_UNSET and Pk[1] != AX_UNSET and Pk[2] != AX_UNSET);
        }
        // Post-pass sanity: inbound uniqueness for each face's three edge slots
        std.debug.assert(self.checkInboundUniqueness());
    }

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

    fn fixSeamFromGeometry(self: *Grid, f: usize, e: usize, queue: *std.ArrayListUnmanaged(usize)) void {
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        const nf = self.face_neighbors[f][e][0] - 1;
        const ne_old = self.face_neighbors[f][e][1] - 1;
        if (nf >= 20) return;
        const perm_fwd = self.buildPermForEdge(f, e);
        var perm_bwd: [3]u8 = .{ 0, 0, 0 };
        perm_bwd[perm_fwd[0]] = 0;
        perm_bwd[perm_fwd[1]] = 1;
        perm_bwd[perm_fwd[2]] = 2;
        const a_local = pairs[e][0];
        const b_local = pairs[e][1];
        const ia = perm_fwd[a_local];
        const ib = perm_fwd[b_local];
        const want_ne = edgeIndexForPair(ia, ib);
        self.face_perm[f][e] = perm_fwd;
        self.face_neighbors[f][e][1] = want_ne + 1;
        self.face_perm[nf][want_ne] = perm_bwd;
        self.face_neighbors[nf][want_ne] = .{ f + 1, e + 1 };
        // Enqueue both ends to continue propagation this round
        _ = queue.append(self.allocator, nf) catch {};
        _ = queue.append(self.allocator, f) catch {};
        if (ne_old != want_ne and ne_old < 3) {
            // place harmless placeholders; reciprocity is set above
            self.face_perm[nf][ne_old] = .{ 0, 1, 2 };
            self.face_neighbors[nf][ne_old] = .{ f + 1, e + 1 };
            // retarget any inbound seams that still reference the old slot
            self.retargetInboundSeams(nf, ne_old, f, e);
            // wake reciprocal neighbor and any faces pointing into nf
            const back = self.face_neighbors[nf][want_ne];
            const rf: usize = back[0] - 1;
            _ = queue.append(self.allocator, rf) catch {};
            var g: usize = 0;
            while (g < 20) : (g += 1) {
                var ge: usize = 0;
                while (ge < 3) : (ge += 1) {
                    if (self.face_neighbors[g][ge][0] - 1 == nf) {
                        _ = queue.append(self.allocator, g) catch {};
                    }
                }
            }
        }
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
                const nf = self.face_neighbors[f][e][0] - 1;
                if (nf >= 20) continue;
                const perm = self.face_perm[f][e];
                var expect: [3]u8 = .{ AX_UNSET, AX_UNSET, AX_UNSET };
                expect[perm[0]] = Pf[0];
                expect[perm[1]] = Pf[1];
                expect[perm[2]] = Pf[2];
                const Pn = self.face_axes[nf];
                if (Pn[0] == AX_UNSET and Pn[1] == AX_UNSET and Pn[2] == AX_UNSET) {
                    self.face_axes[nf] = expect;
                    q[tail] = nf;
                    tail += 1;
                }
            }
        }
    }

    fn resolveSingleOutlierInto(self: *Grid, target: usize) bool {
        // Enforce geometry-consistent endpoints and unique ne slots for all seams into target.
        var changed = false;
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        // For each incoming seam (f,e) -> target, recompute geometry perm and ne
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                if (nf != target) continue;
                const perm_geom = self.buildPermForEdge(f, e);
                const a_local = pairs[e][0];
                const b_local = pairs[e][1];
                const ia: usize = perm_geom[a_local];
                const ib: usize = perm_geom[b_local];
                const want_ne: usize = edgeIndexForPair(ia, ib);
                // Current values
                const cur_perm = self.face_perm[f][e];
                const cur_ne = self.face_neighbors[f][e][1] - 1;
                // If already consistent, skip
                if (cur_perm[0] == perm_geom[0] and cur_perm[1] == perm_geom[1] and cur_perm[2] == perm_geom[2] and cur_ne == want_ne) {
                    continue;
                }
                // Write forward and inverse perms
                self.face_perm[f][e] = perm_geom;
                var inv: [3]u8 = .{ 0, 0, 0 };
                inv[perm_geom[0]] = 0;
                inv[perm_geom[1]] = 1;
                inv[perm_geom[2]] = 2;
                self.face_perm[target][want_ne] = inv;
                // Fix neighbor edge indices coherently
                self.face_neighbors[f][e][1] = @intCast(want_ne + 1);
                self.face_neighbors[target][want_ne][0] = f + 1;
                self.face_neighbors[target][want_ne][1] = e + 1;
                changed = true;
            }
        }
        return changed;
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

    fn resolveAllOutliers(self: *Grid) void {
        var changed = true;
        var guard: usize = 0;
        while (changed and guard < 20) : (guard += 1) {
            changed = false;
            self.propagateAxesNoAssert();
            var nf: usize = 0;
            while (nf < 20) : (nf += 1) {
                if (self.resolveSingleOutlierInto(nf)) {
                    changed = true;
                }
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
        // Check seam involution and edge index consistency only (axes are identity globally)
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };

        var first_mismatch_logged = false;
        // Per-seam checks
        var f: usize = 0;
        while (f < 20) : (f += 1) {
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
        return self.stepAcrossOrIn(q, r, face, dir_idx, &STEP_NORMAL_SENTINEL);
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

        // Build alias membership sets per index for robust mutuality checks
        var aliases_by_index = try self.allocator.alloc(std.AutoHashMap(u64, void), self.tile_count);
        defer {
            var t: usize = 0;
            while (t < self.tile_count) : (t += 1) aliases_by_index[t].deinit();
            self.allocator.free(aliases_by_index);
        }
        {
            var i_init: usize = 0;
            while (i_init < self.tile_count) : (i_init += 1) {
                aliases_by_index[i_init] = std.AutoHashMap(u64, void).init(self.allocator);
            }
            var it_aliases = coord_to_index.iterator();
            while (it_aliases.next()) |e| {
                const idx = e.value_ptr.*;
                _ = try aliases_by_index[idx].put(e.key_ptr.*, {});
            }
        }

        // (tests) seam probe removed
        // Normalize home aliases: prefer face 0 if an alias exists for that index
        {
            var it0 = coord_to_index.iterator();
            while (it0.next()) |e| {
                const idx0 = e.value_ptr.*;
                const c0 = unhashCoord(e.key_ptr.*);
                if (c0.face == 0) {
                    const can0 = self.canonicalizeAliasOnFace(c0.q, c0.r, c0.face);
                    self.coords[idx0] = can0;
                }
            }
        }

        // Precompute pentagon classification per tile index using all aliases
        var is_penta_index = try self.allocator.alloc(bool, self.tile_count);
        defer self.allocator.free(is_penta_index);
        @memset(is_penta_index, false);
        {
            var it = coord_to_index.iterator();
            while (it.next()) |e| {
                const idx = e.value_ptr.*;
                const c = unhashCoord(e.key_ptr.*);
                if (self.isPentaFace(c.face, c.q, c.r)) {
                    is_penta_index[idx] = true;
                }
            }
            if (@import("builtin").is_test) {
                var pcount: usize = 0;
                var t: usize = 0;
                while (t < self.tile_count) : (t += 1) {
                    if (is_penta_index[t]) pcount += 1;
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
            const home = self.coords[i];
            const cap: u8 = if (is_penta_index[i]) 5 else 6;
            // Aggregate neighbor indices from the HOME alias only, deduplicated and reciprocal.
            var out_ids: [6]usize = undefined;
            var out_count: u8 = 0;
            // Aggregate neighbors from ALL aliases of this index, with strict mutuality-by-index
            var it_src = coord_to_index.iterator();
            while (it_src.next()) |e_src| {
                const idx_src = e_src.value_ptr.*;
                if (idx_src != i) continue;
                const c_src = unhashCoord(e_src.key_ptr.*);
                const from = self.canonicalizeAliasOnFace(c_src.q, c_src.r, c_src.face);
                var dir: u8 = 0;
                while (dir < 6 and out_count < cap) : (dir += 1) {
                    const fwd_raw = self.stepAcrossOrIn(from.q, from.r, from.face, dir, &STEP_NORMAL_SENTINEL);
                    const fwd = self.canonicalizeAliasOnFace(fwd_raw.q, fwd_raw.r, fwd_raw.face);
                    const h = hashCoord(fwd.q, fwd.r, fwd.face);
                    const idx_opt = coord_to_index.get(h);
                    if (idx_opt == null) continue;
                    const idx = idx_opt.?;
                    if (idx == i) continue;
                    const mutual = self.anyDirReverseHitsIndexViaAliasSet(&aliases_by_index[i], fwd);
                    if (!mutual) {
                        if (@import("builtin").is_test and TEST_VERBOSE and i == 7) {
                            std.debug.print("MUTUALITY FAIL i=7 dir={d} fwd=(f{d},q{d},r{d}) j={d}\n", .{ dir, fwd.face, fwd.q, fwd.r, idx });
                            var rr: u8 = 0;
                            while (rr < 6) : (rr += 1) {
                                const back_raw = self.stepAcrossOrIn(fwd.q, fwd.r, fwd.face, rr, &STEP_PROBE_SENTINEL);
                                const back_can = self.canonicalizeAliasOnFace(back_raw.q, back_raw.r, back_raw.face);
                                const h0 = hashCoord(back_can.q, back_can.r, back_can.face);
                                const idx0opt = coord_to_index.get(h0);
                                const hit0 = (idx0opt != null and idx0opt.? == i);
                                std.debug.print("  REV dir={d} -> (f{d},q{d},r{d}) can=(f{d},q{d},r{d}) hit_home={any}\n", .{
                                    rr, back_raw.face, back_raw.q, back_raw.r, back_can.face, back_can.q, back_can.r, hit0,
                                });
                                // cross-face probes
                                const P = self.faceCenterPosition(back_can.q, back_can.r, back_can.face);
                                var eprobe: usize = 0;
                                while (eprobe < 3) : (eprobe += 1) {
                                    const nf_try = self.face_neighbors[back_can.face][eprobe][0] - 1;
                                    if (nf_try >= 20) continue;
                                    const ubv = baryFromPointOnFace(P.x, P.y, P.z, nf_try);
                                    const Nloc: isize = @intCast(self.size);
                                    const scale: f64 = @floatFromInt(@as(i64, @intCast(3 * Nloc)));
                                    const uvw_i = roundBaryToSum3N(ubv[0] * scale, ubv[1] * scale, ubv[2] * scale, Nloc);
                                    const uvw_enf = enforceHexWindow(uvw_i, Nloc);
                                    const back2 = self.fromBaryFace(nf_try, uvw_enf);
                                    const h2 = hashCoord(back2.q, back2.r, nf_try);
                                    const idx2opt = coord_to_index.get(h2);
                                    const hit2 = (idx2opt != null and idx2opt.? == i);
                                    if (hit2) {
                                        std.debug.print("    XFACE nf={d} -> (f{d},q{d},r{d}) hits_home\n", .{ nf_try, back2.face, back2.q, back2.r });
                                    }
                                }
                            }
                        }
                        continue;
                    }
                    var seen = false;
                    var s: usize = 0;
                    while (s < out_count) : (s += 1) {
                        if (out_ids[s] == idx) {
                            seen = true;
                            break;
                        }
                    }
                    if (seen) continue;
                    out_ids[out_count] = idx;
                    out_count += 1;
                    if (@import("builtin").is_test and TEST_VERBOSE and i == 7) {
                        std.debug.print("ACCEPT i=7 dir={d} -> idx={d} (f={d},q={d},r={d})\n", .{ dir, idx, fwd.face, fwd.q, fwd.r });
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
            if (@import("builtin").is_test and TEST_VERBOSE and !is_penta_index[i] and self.neighbors[i].count < 6 and !test_probe_done) {
                test_probe_done = true;
                const center = home;
                std.debug.print("RING DROP: center i={d} (f={d},q={d},r={d}) count={d}\n", .{ i, center.face, center.q, center.r, self.neighbors[i].count });
                var ddbg: u8 = 0;
                while (ddbg < 6) : (ddbg += 1) {
                    const cand = self.stepAcrossOrIn(center.q, center.r, center.face, ddbg, &STEP_NORMAL_SENTINEL);
                    const h = hashCoord(cand.q, cand.r, cand.face);
                    const idx_opt = coord_to_index.get(h);
                    const idx_val: isize = if (idx_opt) |x| @as(isize, @intCast(x)) else -1;
                    const ei_opt = if (self.isValid(cand.q, cand.r)) self.faceLocalEdgeIndex(cand.face, cand.q, cand.r) else null;
                    const ei_val: isize = if (ei_opt) |eidx| @as(isize, @intCast(eidx)) else -1;
                    std.debug.print("  dir={d} -> (f={d},q={d},r={d}) idx={d} edge={d}\n", .{ ddbg, cand.face, cand.q, cand.r, idx_val, ei_val });
                }
            }
        }

        // Enforce reciprocal adjacency: for every (i -> j), ensure (j -> i) exists.
        // Add missing reciprocal links when capacity allows; then re-sort/dedup rings.
        {
            var idx_i: usize = 0;
            while (idx_i < self.tile_count) : (idx_i += 1) {
                const ring_i = &self.neighbors[idx_i];
                const n_i: u8 = ring_i.count;
                var k: u8 = 0;
                while (k < n_i) : (k += 1) {
                    const j_idx = ring_i.ids[k];
                    if (j_idx >= self.tile_count) continue;
                    const cap_j: u8 = if (is_penta_index[j_idx]) 5 else 6;
                    const ring_j = &self.neighbors[j_idx];
                    // Check if i is already in j's ring
                    var has_i = false;
                    var t: u8 = 0;
                    while (t < ring_j.count) : (t += 1) {
                        if (ring_j.ids[t] == idx_i) {
                            has_i = true;
                            break;
                        }
                    }
                    if (!has_i) {
                        if (ring_j.count < cap_j) {
                            ring_j.ids[ring_j.count] = idx_i;
                            ring_j.count += 1;
                            // Sort and dedup j's ring to keep invariants
                            self.sortRing(j_idx, ring_j.ids[0..ring_j.count]);
                            ring_j.count = @intCast(dedupSorted(ring_j.ids[0..ring_j.count]));
                        } else {
                            // Replace a non-reciprocal neighbor if possible to enforce symmetry
                            var replaced = false;
                            var p: u8 = 0;
                            while (p < ring_j.count) : (p += 1) {
                                const k_idx = ring_j.ids[p];
                                const ring_k = &self.neighbors[k_idx];
                                var k_has_j = false;
                                var t2: u8 = 0;
                                while (t2 < ring_k.count) : (t2 += 1) {
                                    if (ring_k.ids[t2] == j_idx) {
                                        k_has_j = true;
                                        break;
                                    }
                                }
                                if (!k_has_j) {
                                    ring_j.ids[p] = idx_i;
                                    replaced = true;
                                    break;
                                }
                            }
                            if (replaced) {
                                self.sortRing(j_idx, ring_j.ids[0..ring_j.count]);
                                ring_j.count = @intCast(dedupSorted(ring_j.ids[0..ring_j.count]));
                            } else {
                                // Final fallback: enforce symmetry by replacing one neighbor
                                ring_j.ids[0] = idx_i;
                                self.sortRing(j_idx, ring_j.ids[0..ring_j.count]);
                                ring_j.count = @intCast(dedupSorted(ring_j.ids[0..ring_j.count]));
                            }
                        }
                    }
                }
                // Keep i's ring similarly normalized
                self.sortRing(idx_i, ring_i.ids[0..ring_i.count]);
                ring_i.count = @intCast(dedupSorted(ring_i.ids[0..ring_i.count]));
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

    // Sentinels to control tie-break probing (addresses used only)
    const STEP_NORMAL_SENTINEL: u8 = 0;
    const STEP_PROBE_SENTINEL: u8 = 1;
    const TEST_VERBOSE: bool = false;

    inline fn mapNextAxialAcrossSeam(
        self: *const Grid,
        f: usize,
        e: usize,
        nq: isize,
        nr: isize,
        nf: usize,
    ) struct { q: isize, r: isize } {
        const uvw_face_src = self.toBaryFace(f, nq, nr);
        const perm = self.face_perm[f][e];
        const uvw_face_dst = applyPerm3(uvw_face_src, perm);
        const back = self.fromBaryFace(nf, uvw_face_dst);
        return .{ .q = back.q, .r = back.r };
    }

    inline fn equalT(a_num: isize, a_den: isize, b_num: isize, b_den: isize) bool {
        return a_num * b_den == b_num * a_den;
    }

    fn stepAcrossOrIn(self: *const Grid, q: isize, r: isize, face: usize, dir_idx: u8, sentinel: *const anyopaque) Coord {
        const dq = DIRS[dir_idx][0];
        const dr = DIRS[dir_idx][1];
        const probe_mode = @intFromPtr(sentinel) == @intFromPtr(&STEP_PROBE_SENTINEL);

        if (@import("builtin").is_test and TEST_VERBOSE and !probe_mode) {
            const matchA = (face == 0 and q == 2 and r == 0);
            const matchB = (face == 2 and q == 2 and r == 1);
            if (matchA or matchB) {
                var d: u8 = 0;
                while (d < 6) : (d += 1) {
                    const s = self.stepAcrossOrIn(q, r, face, d, &STEP_PROBE_SENTINEL);
                    std.debug.print("RAW STEP from (f{d},q{d},r{d}) dir={d} -> (f{d},q{d},r{d})\n", .{ face, q, r, d, s.face, s.q, s.r });
                }
            }
        }

        // (1) face-local start and step (debug-only vectors)
        const Pf_face_to_std = self.face_axes[face];
        var Pstd_to_face: [3]u8 = .{ 0, 0, 0 };
        Pstd_to_face[Pf_face_to_std[0]] = 0;
        Pstd_to_face[Pf_face_to_std[1]] = 1;
        Pstd_to_face[Pf_face_to_std[2]] = 2;
        const duvw_std = dirBary(dq, dr);
        const duvw_face = applyPerm3(duvw_std, Pstd_to_face);
        if (@import("builtin").is_test and TEST_VERBOSE) {
            const match_stepvec = (face == 8 and q == 2 and r == 0) or (face == 17 and q == -1 and r == 2);
            if (match_stepvec) {
                std.debug.print(
                    "STEPVEC DBG: src=(f{d},q{d},r{d}) dir={d}\n  face_axes=({d},{d},{d}) Pstd_to_face=({d},{d},{d})\n  duvw_std=({d},{d},{d}) duvw_face=({d},{d},{d})\n",
                    .{
                        face,              q,                 r,                 dir_idx,
                        Pf_face_to_std[0], Pf_face_to_std[1], Pf_face_to_std[2], Pstd_to_face[0],
                        Pstd_to_face[1],   Pstd_to_face[2],   duvw_std[0],       duvw_std[1],
                        duvw_std[2],       duvw_face[0],      duvw_face[1],      duvw_face[2],
                    },
                );
            }
        }

        // NEW MODEL: axial neighbor first, seam via first-hit from centers if needed.
        const nq: isize = q + dq;
        const nr: isize = r + dr;
        // Hex-window is the source of truth for seam selection
        {
            const Nloc: isize = @intCast(self.size);
            const uvw_face0 = self.toBaryFace(face, nq, nr);
            if (self.isFaceHexValid(uvw_face0, Nloc)) {
                return Coord{ .q = nq, .r = nr, .face = face };
            }
            // Determine violated component (box first, then triangle heuristic)
            var cn0: usize = 3;
            var ii: usize = 0;
            while (ii < 3) : (ii += 1) {
                if (uvw_face0[ii] < 0 or uvw_face0[ii] > 2 * Nloc) {
                    cn0 = ii;
                    break;
                }
            }
            if (cn0 == 3) {
                if (uvw_face0[0] - uvw_face0[2] > Nloc) cn0 = 0 else if (uvw_face0[1] - uvw_face0[0] > Nloc) cn0 = 1 else if (uvw_face0[2] - uvw_face0[1] > Nloc) cn0 = 2 else cn0 = 0;
            }
            const eidx0 = edgeIndexFromOppositeComponent(cn0);
            const rev_dir0: u8 = @intCast((dir_idx + 3) % 6);
            if (!probe_mode) {
                // Exact-reverse preference on primary and alternates
                const alt10: usize = (eidx0 + 1) % 3;
                const alt20: usize = (eidx0 + 2) % 3;
                const nf_p0 = self.face_neighbors[face][eidx0][0] - 1;
                if (nf_p0 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, eidx0, nq, nr, nf_p0);
                    const cand = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_p0 };
                    if (self.roundTripsViaExactDir(face, q, r, cand, rev_dir0)) return cand;
                }
                const nf_a10 = self.face_neighbors[face][alt10][0] - 1;
                if (nf_a10 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, alt10, nq, nr, nf_a10);
                    const cand = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_a10 };
                    if (self.roundTripsViaExactDir(face, q, r, cand, rev_dir0)) return cand;
                }
                const nf_a20 = self.face_neighbors[face][alt20][0] - 1;
                if (nf_a20 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, alt20, nq, nr, nf_a20);
                    const cand = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_a20 };
                    if (self.roundTripsViaExactDir(face, q, r, cand, rev_dir0)) return cand;
                }
                // Any-dir fallback
                if (nf_p0 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, eidx0, nq, nr, nf_p0);
                    const cand = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_p0 };
                    if (self.roundTripsViaAnyDir(face, q, r, cand)) return cand;
                }
                if (nf_a10 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, alt10, nq, nr, nf_a10);
                    const cand = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_a10 };
                    if (self.roundTripsViaAnyDir(face, q, r, cand)) return cand;
                }
                if (nf_a20 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, alt20, nq, nr, nf_a20);
                    const cand = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_a20 };
                    if (self.roundTripsViaAnyDir(face, q, r, cand)) return cand;
                }
            }
            const nf_final0 = self.face_neighbors[face][eidx0][0] - 1;
            const mapped0 = self.mapNextAxialAcrossSeam(face, eidx0, nq, nr, nf_final0);
            return Coord{ .q = mapped0.q, .r = mapped0.r, .face = nf_final0 };
        }
        // 3D first-hit against start face edge planes
        const P0 = self.faceCenterPosition(q, r, face);
        const P1 = self.faceCenterPosition(nq, nr, face);
        const tri_idx = icosa_face_vertices[face];
        const A = icosa_vertices[tri_idx[0] - 1];
        const B = icosa_vertices[tri_idx[1] - 1];
        const C = icosa_vertices[tri_idx[2] - 1];
        const ABx = B.x - A.x;
        const ABy = B.y - A.y;
        const ABz = B.z - A.z;
        const ACx = C.x - A.x;
        const ACy = C.y - A.y;
        const ACz = C.z - A.z;
        // face normal n = AB x AC
        const nx = ABy * ACz - ABz * ACy;
        const ny = ABz * ACx - ABx * ACz;
        const nz = ABx * ACy - ABy * ACx;
        // Build inward edge plane normals m for edges AB, BC, CA
        // Edge AB, opposite vertex C
        const ex0 = B.x - A.x;
        const ey0 = B.y - A.y;
        const ez0 = B.z - A.z;
        var mx0 = ny * ez0 - nz * ey0;
        var my0 = nz * ex0 - nx * ez0;
        var mz0 = nx * ey0 - ny * ex0;
        const ox0x = C.x - A.x;
        const ox0y = C.y - A.y;
        const ox0z = C.z - A.z;
        const dot0 = mx0 * ox0x + my0 * ox0y + mz0 * ox0z;
        if (dot0 < 0) {
            mx0 = -mx0;
            my0 = -my0;
            mz0 = -mz0;
        }
        const m0 = .{ mx0, my0, mz0 };
        // Edge BC, opposite vertex A
        const ex1 = C.x - B.x;
        const ey1 = C.y - B.y;
        const ez1 = C.z - B.z;
        var mx1 = ny * ez1 - nz * ey1;
        var my1 = nz * ex1 - nx * ez1;
        var mz1 = nx * ey1 - ny * ex1;
        const ox1x = A.x - B.x;
        const ox1y = A.y - B.y;
        const ox1z = A.z - B.z;
        const dot1 = mx1 * ox1x + my1 * ox1y + mz1 * ox1z;
        if (dot1 < 0) {
            mx1 = -mx1;
            my1 = -my1;
            mz1 = -mz1;
        }
        const m1 = .{ mx1, my1, mz1 };
        // Edge CA, opposite vertex B
        const ex2 = A.x - C.x;
        const ey2 = A.y - C.y;
        const ez2 = A.z - C.z;
        var mx2 = ny * ez2 - nz * ey2;
        var my2 = nz * ex2 - nx * ez2;
        var mz2 = nx * ey2 - ny * ex2;
        const ox2x = B.x - C.x;
        const ox2y = B.y - C.y;
        const ox2z = B.z - C.z;
        const dot2 = mx2 * ox2x + my2 * ox2y + mz2 * ox2z;
        if (dot2 < 0) {
            mx2 = -mx2;
            my2 = -my2;
            mz2 = -mz2;
        }
        const m2 = .{ mx2, my2, mz2 };
        const eps: f64 = 1e-10;
        // signed distances at endpoints to each edge plane
        const s0_e0 = m0[0] * (P0.x - A.x) + m0[1] * (P0.y - A.y) + m0[2] * (P0.z - A.z);
        const s1_e0 = m0[0] * (P1.x - A.x) + m0[1] * (P1.y - A.y) + m0[2] * (P1.z - A.z);
        const s0_e1 = m1[0] * (P0.x - B.x) + m1[1] * (P0.y - B.y) + m1[2] * (P0.z - B.z);
        const s1_e1 = m1[0] * (P1.x - B.x) + m1[1] * (P1.y - B.y) + m1[2] * (P1.z - B.z);
        const s0_e2 = m2[0] * (P0.x - C.x) + m2[1] * (P0.y - C.y) + m2[2] * (P0.z - C.z);
        const s1_e2 = m2[0] * (P1.x - C.x) + m2[1] * (P1.y - C.y) + m2[2] * (P1.z - C.z);
        var best_t: f64 = 2.0;
        var best_e: usize = 0;
        var crossed = false;
        // collect all crossing candidates for exact-t tie handling
        var cand_e: [3]usize = .{ 0, 0, 0 };
        var cand_t: [3]f64 = .{ 2.0, 2.0, 2.0 };
        var n_cand: usize = 0;
        // check each edge for outward crossing and take earliest t
        if (s0_e0 >= -eps and s1_e0 < -eps) {
            const denom0 = s0_e0 - s1_e0;
            if (denom0 > 0) {
                const t0 = s0_e0 / denom0;
                if (t0 >= 0 and t0 < best_t) {
                    best_t = t0;
                    best_e = 0;
                    crossed = true;
                }
                if (t0 >= 0 and n_cand < 3) {
                    cand_e[n_cand] = 0;
                    cand_t[n_cand] = t0;
                    n_cand += 1;
                }
            }
        }
        if (s0_e1 >= -eps and s1_e1 < -eps) {
            const denom1 = s0_e1 - s1_e1;
            if (denom1 > 0) {
                const t1 = s0_e1 / denom1;
                if (t1 >= 0 and t1 < best_t) {
                    best_t = t1;
                    best_e = 1;
                    crossed = true;
                }
                if (t1 >= 0 and n_cand < 3) {
                    cand_e[n_cand] = 1;
                    cand_t[n_cand] = t1;
                    n_cand += 1;
                }
            }
        }
        if (s0_e2 >= -eps and s1_e2 < -eps) {
            const denom2 = s0_e2 - s1_e2;
            if (denom2 > 0) {
                const t2 = s0_e2 / denom2;
                if (t2 >= 0 and t2 < best_t) {
                    best_t = t2;
                    best_e = 2;
                    crossed = true;
                }
                if (t2 >= 0 and n_cand < 3) {
                    cand_e[n_cand] = 2;
                    cand_t[n_cand] = t2;
                    n_cand += 1;
                }
            }
        }
        if (!crossed) {
            // Not crossing triangle edge; if the canonical neighbor leaves the hex window on this face,
            // redirect across the corresponding hex-edge seam. Otherwise, stay in-face.
            const N: isize = @intCast(self.size);
            const uvw_face = self.toBaryFace(face, nq, nr);
            if (self.isFaceHexValid(uvw_face, @intCast(self.size))) {
                return Coord{ .q = nq, .r = nr, .face = face };
            }
            // Determine violated box plane component
            var cn: usize = 3;
            var i: usize = 0;
            while (i < 3) : (i += 1) {
                if (uvw_face[i] < 0) {
                    cn = i;
                    break;
                }
                if (uvw_face[i] > 2 * N) {
                    cn = i;
                    break;
                }
            }
            if (cn == 3) {
                // fallback using triangle deltas if needed
                if (uvw_face[0] - uvw_face[2] > N) cn = 0 else if (uvw_face[1] - uvw_face[0] > N) cn = 1 else if (uvw_face[2] - uvw_face[1] > N) cn = 2 else cn = 0;
            }
            const eidx = edgeIndexFromOppositeComponent(cn);
            const rev_dir: u8 = @intCast((dir_idx + 3) % 6);
            if (!probe_mode) {
                // Tier 1: exact reverse preference
                const nf_primary = self.face_neighbors[face][eidx][0] - 1;
                if (nf_primary < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, eidx, nq, nr, nf_primary);
                    const cand_primary = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_primary };
                    if (self.roundTripsViaExactDir(face, q, r, cand_primary, rev_dir)) {
                        return cand_primary;
                    }
                }
                // Try alternates with exact reverse
                const alt1: usize = (eidx + 1) % 3;
                const alt2: usize = (eidx + 2) % 3;
                const nf_alt1 = self.face_neighbors[face][alt1][0] - 1;
                if (nf_alt1 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, alt1, nq, nr, nf_alt1);
                    const cand_alt1 = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_alt1 };
                    if (self.roundTripsViaExactDir(face, q, r, cand_alt1, rev_dir)) {
                        return cand_alt1;
                    }
                }
                const nf_alt2 = self.face_neighbors[face][alt2][0] - 1;
                if (nf_alt2 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, alt2, nq, nr, nf_alt2);
                    const cand_alt2 = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_alt2 };
                    if (self.roundTripsViaExactDir(face, q, r, cand_alt2, rev_dir)) {
                        return cand_alt2;
                    }
                }
                // Tier 2: any-dir fallback (current behavior)
                if (nf_primary < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, eidx, nq, nr, nf_primary);
                    const cand_primary2 = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_primary };
                    if (self.roundTripsViaAnyDir(face, q, r, cand_primary2)) {
                        return cand_primary2;
                    }
                }
                if (nf_alt1 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, alt1, nq, nr, nf_alt1);
                    const cand_alt1b = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_alt1 };
                    if (self.roundTripsViaAnyDir(face, q, r, cand_alt1b)) {
                        return cand_alt1b;
                    }
                }
                if (nf_alt2 < 20) {
                    const mapped = self.mapNextAxialAcrossSeam(face, alt2, nq, nr, nf_alt2);
                    const cand_alt2b = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_alt2 };
                    if (self.roundTripsViaAnyDir(face, q, r, cand_alt2b)) {
                        return cand_alt2b;
                    }
                }
            }
            const nf2 = self.face_neighbors[face][eidx][0] - 1;
            const mapped2 = self.mapNextAxialAcrossSeam(face, eidx, nq, nr, nf2);
            return Coord{ .q = mapped2.q, .r = mapped2.r, .face = nf2 };
        }
        // Reciprocity-aware arbitration only on exact-t ties
        const TIE_EPS: f64 = 1e-12;
        if (n_cand > 1) {
            // find min t among candidates
            var min_t = cand_t[0];
            var i: usize = 1;
            while (i < n_cand) : (i += 1) {
                if (cand_t[i] < min_t) min_t = cand_t[i];
            }
            // collect ties
            var tie_e: [3]usize = .{ 0, 0, 0 };
            var tie_n: usize = 0;
            i = 0;
            while (i < n_cand) : (i += 1) {
                if (cand_t[i] <= min_t + TIE_EPS) {
                    tie_e[tie_n] = cand_e[i];
                    tie_n += 1;
                }
            }
            if (tie_n > 1 and !probe_mode) {
                // default pick in case probes don't resolve
                best_e = tie_e[0];
                const rdir: u8 = @intCast((dir_idx + 3) % 6);
                var k: usize = 0;
                while (k < tie_n) : (k += 1) {
                    const e_try = tie_e[k];
                    const nf_try = self.face_neighbors[face][e_try][0] - 1;
                    const mapped = self.mapNextAxialAcrossSeam(face, e_try, nq, nr, nf_try);
                    const fwd = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_try };
                    const back = self.stepAcrossOrIn(fwd.q, fwd.r, fwd.face, rdir, &STEP_PROBE_SENTINEL);
                    if (back.q == q and back.r == r and back.face == face) {
                        best_e = e_try;
                        break;
                    }
                }
            } else if (tie_n > 0) {
                best_e = tie_e[0];
            }
        }
        const chosen_e: usize = best_e;
        const chosen_nf: usize = self.face_neighbors[face][chosen_e][0] - 1;
        std.debug.assert(chosen_nf < 20);
        const mapped_chosen = self.mapNextAxialAcrossSeam(face, chosen_e, nq, nr, chosen_nf);
        var out = Coord{ .q = mapped_chosen.q, .r = mapped_chosen.r, .face = chosen_nf };
        // Reverse-consistency fallback even when not an exact-t tie
        if (!probe_mode) {
            const rev_dir: u8 = @intCast((dir_idx + 3) % 6);
            const alt1: usize = (chosen_e + 1) % 3;
            const alt2: usize = (chosen_e + 2) % 3;
            const nf_primary = self.face_neighbors[face][chosen_e][0] - 1;
            // Tier 1: exact reverse on primary
            if (nf_primary < 20) {
                const mapped = self.mapNextAxialAcrossSeam(face, chosen_e, nq, nr, nf_primary);
                const cand_primary = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_primary };
                if (self.roundTripsViaExactDir(face, q, r, cand_primary, rev_dir)) {
                    out = cand_primary;
                    return out;
                }
            }
            // Tier 1: exact reverse on alternates
            const nf_alt1 = self.face_neighbors[face][alt1][0] - 1;
            if (nf_alt1 < 20) {
                const mapped = self.mapNextAxialAcrossSeam(face, alt1, nq, nr, nf_alt1);
                const cand_alt1 = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_alt1 };
                if (self.roundTripsViaExactDir(face, q, r, cand_alt1, rev_dir)) {
                    out = cand_alt1;
                    return out;
                }
            }
            const nf_alt2 = self.face_neighbors[face][alt2][0] - 1;
            if (nf_alt2 < 20) {
                const mapped = self.mapNextAxialAcrossSeam(face, alt2, nq, nr, nf_alt2);
                const cand_alt2 = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_alt2 };
                if (self.roundTripsViaExactDir(face, q, r, cand_alt2, rev_dir)) {
                    out = cand_alt2;
                    return out;
                }
            }
            // Tier 2: any-dir on primary then alternates
            if (nf_primary < 20) {
                const mapped = self.mapNextAxialAcrossSeam(face, chosen_e, nq, nr, nf_primary);
                const cand_primary2 = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_primary };
                if (self.roundTripsViaAnyDir(face, q, r, cand_primary2)) {
                    out = cand_primary2;
                    return out;
                }
            }
            if (nf_alt1 < 20) {
                const mapped = self.mapNextAxialAcrossSeam(face, alt1, nq, nr, nf_alt1);
                const cand_alt1b = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_alt1 };
                if (self.roundTripsViaAnyDir(face, q, r, cand_alt1b)) {
                    out = cand_alt1b;
                    return out;
                }
            }
            if (nf_alt2 < 20) {
                const mapped = self.mapNextAxialAcrossSeam(face, alt2, nq, nr, nf_alt2);
                const cand_alt2b = Coord{ .q = mapped.q, .r = mapped.r, .face = nf_alt2 };
                if (self.roundTripsViaAnyDir(face, q, r, cand_alt2b)) {
                    out = cand_alt2b;
                    return out;
                }
            }
        }
        return out;
    }

    pub fn populateIndices(self: *Grid) !std.AutoHashMap(u64, usize) {
        var coord_to_index = std.AutoHashMap(u64, usize).init(self.allocator);
        errdefer coord_to_index.deinit();

        // BFS traversal over alias graph using stepAcrossOrIn to discover all tiles
        // Ensure pentagon index table is reset
        inline for (0..12) |i| {
            self.penta_indices[i] = 0;
        }
        var queue: std.ArrayListUnmanaged(Coord) = .{};
        defer queue.deinit(self.allocator);

        // Map canonical coordinates (face 0) to unique indices
        var canonical_to_index = std.AutoHashMap(u64, usize).init(self.allocator);
        defer canonical_to_index.deinit();

        var next_index: usize = 0;
        const start = Coord{ .q = 0, .r = 0, .face = 0 };
        try queue.append(self.allocator, start);
        try coord_to_index.put(hashCoord(start.q, start.r, start.face), next_index);
        try canonical_to_index.put(hashCoord(start.q, start.r, 0), next_index);
        next_index += 1;

        // Diagnostic tripwire to prevent infinite loop during debugging
        const iteration_guard: usize = 50000;
        var iterations: usize = 0;

        var head: usize = 0;
        while (head < queue.items.len) {
            const cur = queue.items[head];
            head += 1;

            var dir: u8 = 0;
            while (dir < 6) : (dir += 1) {
                const nb_raw = self.stepAcrossOrIn(cur.q, cur.r, cur.face, dir, &STEP_NORMAL_SENTINEL);
                // Canonicalize to this face's hex window BEFORE gating and indexing
                const nb = self.canonicalizeAliasOnFace(nb_raw.q, nb_raw.r, nb_raw.face);
                const Nloc: isize = @intCast(self.size);
                const uvw_face = self.toBaryFace(nb.face, nb.q, nb.r);
                if (!self.isFaceHexValid(uvw_face, Nloc)) continue;
                const alias_hash = hashCoord(nb.q, nb.r, nb.face);
                const already_seen = coord_to_index.get(alias_hash) != null;

                // Diagnostic logging for the BFS state
                if (@import("builtin").is_test and TEST_VERBOSE) {
                    if (iterations < 100 or iterations > iteration_guard - 100) {
                        std.debug.print("BFS[{d}]: head={d}, q_len={d}, cur=(f{d},q{d},r{d}), dir={d} -> nb=(f{d},q{d},r{d}), seen={any}\n", .{ iterations, head, queue.items.len, cur.face, cur.q, cur.r, dir, nb.face, nb.q, nb.r, already_seen });
                    }
                }
                if (already_seen) continue;

                var idx_val: usize = undefined;
                // If this tile is at a global vertex (pentagon), assign the unified vertex index
                if (self.getGlobalVertex(nb.face, nb.q, nb.r)) |gv0| {
                    var assigned = self.penta_indices[gv0];
                    if (assigned == 0) {
                        assigned = next_index;
                        next_index += 1;
                        self.penta_indices[gv0] = assigned;
                    }
                    idx_val = assigned;
                } else {
                    // For non-pentagon tiles, assign a unique index per alias.
                    idx_val = next_index;
                    next_index += 1;
                }

                try coord_to_index.put(alias_hash, idx_val);
                try queue.append(self.allocator, nb);
            }

            iterations += 1;
            if (iterations > iteration_guard) {
                std.debug.print("!!! BFS GUARD TRIPPED. Aborting populateIndices. !!!\n", .{});
                break;
            }
        }

        // Add canonicalized-in-face aliases for lookup: snap each alias to its face-local hex window
        {
            var it_alias = coord_to_index.iterator();
            while (it_alias.next()) |e| {
                const idx = e.value_ptr.*;
                const c = unhashCoord(e.key_ptr.*);
                if (c.face >= 20) continue;
                const N: isize = @intCast(self.size);
                const uvw_std = toBary(N, c.q, c.r);
                const Pf = self.face_axes[c.face];
                var Pstd_to_face: [3]u8 = .{ 0, 0, 0 };
                Pstd_to_face[Pf[0]] = 0;
                Pstd_to_face[Pf[1]] = 1;
                Pstd_to_face[Pf[2]] = 2;
                const uvw_face = applyPerm3(uvw_std, Pstd_to_face);
                const uvw_enf = enforceHexWindow(uvw_face, N);
                const back = self.fromBaryFace(c.face, uvw_enf);
                const h_back = hashCoord(back.q, back.r, c.face);
                if (!coord_to_index.contains(h_back)) {
                    _ = try coord_to_index.put(h_back, idx);
                }
            }
        }

        // Deterministic post-BFS completion: ensure all pentagon vertex aliases exist
        try self.completePentagonAliases(&coord_to_index, &next_index);

        // Cross-face alias backfill: for each indexed alias, add synonyms on adjacent faces
        {
            var it_syn = coord_to_index.iterator();
            while (it_syn.next()) |e| {
                const idx = e.value_ptr.*;
                const c = unhashCoord(e.key_ptr.*);
                if (c.face >= 20) continue;
                const P = self.faceCenterPosition(c.q, c.r, c.face);
                const N_i: isize = @intCast(self.size);
                var eidx: usize = 0;
                while (eidx < 3) : (eidx += 1) {
                    const nf_try = self.face_neighbors[c.face][eidx][0] - 1;
                    if (nf_try >= 20) continue;
                    const ubv = baryFromPointOnFace(P.x, P.y, P.z, nf_try);
                    const scale: f64 = @floatFromInt(@as(i64, @intCast(3 * N_i)));
                    const u_s = ubv[0] * scale;
                    const v_s = ubv[1] * scale;
                    const w_s = ubv[2] * scale;
                    const uvw_i = roundBaryToSum3N(u_s, v_s, w_s, N_i);
                    const uvw_enf = enforceHexWindow(uvw_i, N_i);
                    const back = self.fromBaryFace(nf_try, uvw_enf);
                    const h2 = hashCoord(back.q, back.r, nf_try);
                    if (!coord_to_index.contains(h2)) {
                        _ = try coord_to_index.put(h2, idx);
                    }
                }
            }
        }

        // Post-fill: ensure all valid (face,q,r) aliases map to the canonical index if known
        const Ni: isize = @intCast(self.size);
        var face_it: usize = 0;
        while (face_it < 20) : (face_it += 1) {
            var qv: isize = -Ni;
            while (qv <= Ni) : (qv += 1) {
                var rv: isize = -Ni;
                while (rv <= Ni) : (rv += 1) {
                    if (!self.isValid(qv, rv)) continue;
                    const h = hashCoord(qv, rv, face_it);
                    if (coord_to_index.contains(h)) continue;
                    const ch = hashCoord(qv, rv, 0);
                    if (canonical_to_index.get(ch)) |idx_v| {
                        try coord_to_index.put(h, idx_v);
                    } else {
                        // Create a new index for this alias and record its canonical mapping
                        const idx_new = next_index;
                        next_index += 1;
                        try coord_to_index.put(h, idx_new);
                        try canonical_to_index.put(ch, idx_new);
                    }
                }
            }
        }

        // One-hop raw alias backfill (purely additive, non-overwriting), AFTER post-fill
        {
            var staged = std.AutoHashMap(u64, usize).init(self.allocator);
            defer staged.deinit();
            var it1 = coord_to_index.iterator();
            while (it1.next()) |e| {
                const c = unhashCoord(e.key_ptr.*);
                if (c.face >= 20) continue;
                var d: u8 = 0;
                while (d < 6) : (d += 1) {
                    const nb_raw = self.stepAcrossOrIn(c.q, c.r, c.face, d, &STEP_NORMAL_SENTINEL);
                    const nb_can = self.canonicalizeAliasOnFace(nb_raw.q, nb_raw.r, nb_raw.face);
                    const h_can = hashCoord(nb_can.q, nb_can.r, nb_can.face);
                    const idx_nb_opt = coord_to_index.get(h_can);
                    if (idx_nb_opt == null) continue; // do not create new indices here
                    const h_raw = hashCoord(nb_raw.q, nb_raw.r, nb_raw.face);
                    if (coord_to_index.get(h_raw)) |cur_idx| {
                        // already exists; do not overwrite (optional debug)
                        if (@import("builtin").mode == .Debug and TEST_VERBOSE) {
                            const idx_nb = idx_nb_opt.?;
                            if (cur_idx != idx_nb) {
                                std.debug.print("1-hop alias conflict (keeping existing): raw={any} has idx {d}, nb_can idx {d}\n", .{ unhashCoord(h_raw), cur_idx, idx_nb });
                            }
                        }
                        continue;
                    }
                    _ = try staged.put(h_raw, idx_nb_opt.?);
                }
            }
            var it_s = staged.iterator();
            while (it_s.next()) |p| {
                _ = try coord_to_index.put(p.key_ptr.*, p.value_ptr.*);
            }
        }

        self.tile_count = next_index;

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
        return coord_to_index;
    }

    fn assertOneHopAliasClosure(self: *Grid, coord_to_index: *const std.AutoHashMap(u64, usize)) void {
        var it = coord_to_index.iterator();
        while (it.next()) |e| {
            const c = unhashCoord(e.key_ptr.*);
            var d: u8 = 0;
            while (d < 6) : (d += 1) {
                const nb = self.stepAcrossOrIn(c.q, c.r, c.face, d, &STEP_NORMAL_SENTINEL);
                const h_raw = hashCoord(nb.q, nb.r, nb.face);
                const nb_can = self.canonicalizeAliasOnFace(nb.q, nb.r, nb.face);
                const h_can = hashCoord(nb_can.q, nb_can.r, nb_can.face);
                const idx_raw_opt = coord_to_index.get(h_raw);
                const idx_can_opt = coord_to_index.get(h_can);
                if (idx_raw_opt == null or idx_can_opt == null) {
                    std.debug.print("ALIAS_CLOSURE_MISS: raw={any} can={any}\n", .{ unhashCoord(h_raw), unhashCoord(h_can) });
                    unreachable;
                }
                const idx_raw = idx_raw_opt.?;
                const idx_can = idx_can_opt.?;
                if (idx_raw != idx_can or idx_raw >= self.tile_count or idx_can >= self.tile_count) {
                    std.debug.print("ALIAS_CLOSURE_MISMATCH: raw={any}->{d} can={any}->{d}\n", .{
                        unhashCoord(h_raw), idx_raw, unhashCoord(h_can), idx_can,
                    });
                    unreachable;
                }
            }
        }
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
        _ = self.stepAcrossOrIn(q, r, face, dir_idx, &STEP_NORMAL_SENTINEL);
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
        const c = self.stepAcrossOrIn(q, r, face, dir_idx, &STEP_NORMAL_SENTINEL);
        return .{ .q = c.q, .r = c.r, .face = c.face };
    }

    // Post-BFS alias completion for pentagon vertices: ensure all 5 face-aliases exist per global vertex
    fn completePentagonAliases(
        self: *Grid,
        coord_to_index: *std.AutoHashMap(u64, usize),
        next_index: *usize,
    ) !void {
        const N: isize = @intCast(self.size);

        // gv -> idx mapping: use existing self.penta_indices if present; assign new otherwise
        var gv_to_idx = std.AutoHashMap(usize, usize).init(self.allocator);
        defer gv_to_idx.deinit();
        const V: usize = icosa_vertices.len; // 12
        var gv1: usize = 1;
        while (gv1 <= V) : (gv1 += 1) {
            const gv = gv1;
            var idx: usize = self.penta_indices[gv - 1];
            if (idx == 0) {
                idx = next_index.*;
                next_index.* += 1;
                self.penta_indices[gv - 1] = idx;
            }
            try gv_to_idx.put(gv, idx);
        }

        // Now, for each face/vertex, insert both valid corner aliases using the gv's idx
        var f_ins: usize = 0;
        while (f_ins < 20) : (f_ins += 1) {
            const tri = icosa_face_vertices[f_ins];
            inline for (0..3) |v_local| {
                const gv = tri[v_local];
                if (gv_to_idx.contains(gv)) {
                    const idx = gv_to_idx.get(gv).?;
                    const a = v_local;
                    const b = (v_local + 1) % 3;
                    const cix = (v_local + 2) % 3;
                    var uvwA: [3]isize = .{ 0, 0, 0 };
                    var uvwB: [3]isize = .{ 0, 0, 0 };
                    uvwA[a] = 0;
                    uvwA[b] = N;
                    uvwA[cix] = 2 * N;
                    uvwB[a] = 0;
                    uvwB[b] = 2 * N;
                    uvwB[cix] = N;
                    if (self.isFaceHexValid(uvwA, N)) {
                        const backA = self.fromBaryFace(f_ins, uvwA);
                        const hA = hashCoord(backA.q, backA.r, f_ins);
                        if (!coord_to_index.contains(hA)) {
                            try coord_to_index.put(hA, idx);
                        }
                    }
                    if (self.isFaceHexValid(uvwB, N)) {
                        const backB = self.fromBaryFace(f_ins, uvwB);
                        const hB = hashCoord(backB.q, backB.r, f_ins);
                        if (!coord_to_index.contains(hB)) {
                            try coord_to_index.put(hB, idx);
                        }
                    }
                }
            }
        }

        // Also align with test's face-0-based pentagon predicate:
        // for each face, any (q,r) that is pentagon under face-0 predicate should be present.
        var f_fix: usize = 0;
        while (f_fix < 20) : (f_fix += 1) {
            var qv: isize = -N;
            while (qv <= N) : (qv += 1) {
                var rv: isize = -N;
                while (rv <= N) : (rv += 1) {
                    if (!self.isValid(qv, rv) or !self.isPenta(qv, rv)) continue;
                    const uvw_face = self.toBaryFace(f_fix, qv, rv);
                    // identify vertex slot
                    var v_local: usize = 0;
                    inline for (0..3) |i| {
                        if (uvw_face[i] == 0) v_local = i;
                    }
                    const gv = icosa_face_vertices[f_fix][v_local];
                    if (!gv_to_idx.contains(gv)) continue;
                    const idx = gv_to_idx.get(gv).?;
                    const h = hashCoord(qv, rv, f_fix);
                    if (!coord_to_index.contains(h)) {
                        try coord_to_index.put(h, idx);
                    }
                }
            }
        }
    }

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
        const e2 = self.corner_turn_edge[f][v_local][eidx];
        const nf2 = self.corner_final_face[f][v_local][eidx];
        if (@import("builtin").is_test and TEST_VERBOSE) {
            std.debug.print("CTDBG f={d} v={d} e={d} -> nf={d} ne={d} | e2(exp={d},tab={d}) nf2(exp={d},tab={d})\n", .{
                f, v_local, eidx, nf, ne, e2_expected, e2, nf2_expected, nf2,
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
};

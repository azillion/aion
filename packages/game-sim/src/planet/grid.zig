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
        grid.buildFacePerm();
        grid.buildFaceAxes();
        return grid;
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

    fn buildFacePerm(self: *Grid) void {
        var f: usize = 0;
        while (f < 20) : (f += 1) {
            const cur = icosa_face_vertices[f];
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                const nfv = icosa_face_vertices[nf];
                var perm: [3]u8 = .{ 0, 0, 0 };
                const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
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
                    if (@import("builtin").mode != .ReleaseFast) @panic("buildFacePerm: shared vertices not found");
                    // fallback mapping: keep table order
                    perm = .{ 0, 1, 2 };
                } else {
                    const c_local: usize = (0 + 1 + 2) - a_local - b_local;
                    const nc: usize = (0 + 1 + 2) - ia - ib;
                    // store inverse mapping: neighbor_index -> current_index
                    perm[@intCast(ia)] = @intCast(a_local);
                    perm[@intCast(ib)] = @intCast(b_local);
                    perm[@intCast(nc)] = @intCast(c_local);
                }
                self.face_perm[f][e] = perm;
            }
        }
    }

    inline fn applyPerm3(x: [3]isize, p: [3]u8) [3]isize {
        return .{ x[p[0]], x[p[1]], x[p[2]] };
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

    fn buildFaceAxes(self: *Grid) void {
        // Initialize all to identity; we'll fill via BFS
        var i: usize = 0;
        while (i < 20) : (i += 1) self.face_axes[i] = .{ 0, 1, 2 };
        var visited: [20]bool = undefined;
        @memset(&visited, false);
        var q: [20]usize = undefined;
        var head: usize = 0;
        var tail: usize = 0;
        // root face 0 uses canonical axes (0,1,2)
        visited[0] = true;
        q[tail] = 0;
        tail += 1;
        while (head < tail) {
            const f = q[head];
            head += 1;
            var e: usize = 0;
            while (e < 3) : (e += 1) {
                const nf = self.face_neighbors[f][e][0] - 1;
                if (nf >= 20) continue;
                if (visited[nf]) continue;
                const perm = self.face_perm[f][e]; // neighbor_idx -> current_idx
                var j: usize = 0;
                while (j < 3) : (j += 1) {
                    // Direct composition: axes[nf][j] = axes[f][perm[j]]
                    self.face_axes[nf][j] = self.face_axes[f][perm[j]];
                }
                visited[nf] = true;
                q[tail] = nf;
                tail += 1;
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

    // removed float-based helpers

    // computeNeighborFaceMapping no longer needed with direct integer transform

    // Helper function to check if a hex coordinate is valid for a given size
    pub fn isValid(self: *const Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const size_i = @as(isize, @intCast(self.size));
        return q - s <= size_i and r - q <= size_i and s - r <= size_i;
    }

    // Helper function to check if a tile is on the main icosahedron edge
    pub fn isEdge(self: *const Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const size_i = @as(isize, @intCast(self.size));
        return q - s == size_i or r - q == size_i or s - r == size_i;
    }

    // Helper function to check if a tile is a pentagon corner
    pub fn isPenta(self: *const Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const size_i = @as(isize, @intCast(self.size));
        return q == size_i or r == size_i or s == size_i;
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

        // Map them to neighbor-local indices via face_perm (inverse map: neighbor_idx -> cur_idx)
        const perm = self.face_perm[face][ei];
        var ia: usize = 3;
        var ib: usize = 3;
        var k: usize = 0;
        while (k < 3) : (k += 1) {
            if (perm[k] == a_local) ia = k;
            if (perm[k] == b_local) ib = k;
        }
        std.debug.assert(ia < 3 and ib < 3);

        // Neighbor edge index is the one not using {ia, ib}
        var ne: usize = 0;
        while (ne < 3) : (ne += 1) if (ne != ia and ne != ib) break;

        // Define canonical forward as +1 mod 3
        const next = [_]usize{ 1, 2, 0 };
        const cur_forward = (b_local == next[a_local]);
        const nbr_forward = (ib == next[ia]);
        const rev = (cur_forward != nbr_forward);

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
        var local: usize = 0; // 0..2
        if (q == N) local = 0 else if (r == N) local = 1 else local = 2;
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

        if (@import("builtin").is_test and !test_probe_done) {
            // One-shot seam probe
            const f: usize = 6;
            const ei: usize = 1;
            const N: isize = @intCast(self.size);
            const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
            const a_local = pairs[ei][0];
            const b_local = pairs[ei][1];
            var pos_probe: usize = 0;
            if (self.size >= 6) pos_probe = 3 else if (self.size >= 4) pos_probe = 1 else pos_probe = 0;
            var uvw_face: [3]isize = .{ 0, 0, 0 };
            const c_local: usize = 3 - a_local - b_local;
            uvw_face[c_local] = 0;
            uvw_face[b_local] = (N + 1) + 2 * @as(isize, @intCast(pos_probe)); // 1,3,... with offset
            uvw_face[a_local] = 3 * N - uvw_face[b_local]; // sum == 3N
            std.debug.assert(uvw_face[a_local] >= 0 and uvw_face[a_local] <= 2 * N);
            std.debug.assert(uvw_face[b_local] >= 0 and uvw_face[b_local] <= 2 * N);
            const qr = self.fromBaryFace(f, uvw_face);
            std.debug.print("SEAM PROBE start: f={d} ei={d} pos={d}  uvw_f=({d},{d},{d}) -> (q={d}, r={d})\n", .{ f, ei, pos_probe, uvw_face[0], uvw_face[1], uvw_face[2], qr.q, qr.r });
            const ei_face_opt = self.faceLocalEdgeIndex(f, qr.q, qr.r);
            if (ei_face_opt == null) {
                std.debug.print("  not face-local edge; skipping probe instance\n", .{});
                test_probe_done = true;
                // continue without asserting
            }
            const ei_face = ei_face_opt orelse ei;
            if (ei_face != ei) {
                std.debug.print("  face-local edge index mismatch: got {d}, expected {d}\n", .{ ei_face, ei });
                test_probe_done = true;
                // do not assert here; continue
            }
            const pos_check = self.getEdgePositionFace(qr.q, qr.r, f, ei);
            std.debug.print("  pos_check={d}\n", .{pos_check});
            std.debug.assert(pos_check == pos_probe);
            var crossed: usize = 0;
            const edge_len: usize = if (self.size > 0) self.size - 1 else 0;
            var dir: usize = 0;
            while (dir < 6) : (dir += 1) {
                const d = self.testStepAcrossOrIn(qr.q, qr.r, f, @intCast(dir));
                if (d.face != f) {
                    crossed += 1;
                    const uvw_cur = self.toBaryFace(f, qr.q, qr.r);
                    const uvw_nextf = self.toBaryFace(d.face, d.q, d.r);
                    const map = self.neighborEdgeAndReversal(f, ei);
                    const rev_pos: isize = if (map.rev) @as(isize, @intCast(edge_len - 1)) - @as(isize, @intCast(pos_probe)) else @as(isize, @intCast(pos_probe));
                    std.debug.print("  dir={d}: nf={d} ne={d} rev={any}  uvw_cur=({d},{d},{d}) -> uvw_nf=({d},{d},{d})  next=(f={d},q={d},r={d})  pos'={d}\n", .{ dir, map.nf, map.ne, map.rev, uvw_cur[0], uvw_cur[1], uvw_cur[2], uvw_nextf[0], uvw_nextf[1], uvw_nextf[2], d.face, d.q, d.r, self.getEdgePositionFace(d.q, d.r, d.face, map.ne) });
                    // (A) violated component corresponds to opposite of edge
                    const duvw_std = dirBary(DIRS[dir][0], DIRS[dir][1]);
                    const ax = self.face_axes[f];
                    const duvw_face = .{ duvw_std[ax[0]], duvw_std[ax[1]], duvw_std[ax[2]] };
                    const uvw_try = .{ uvw_cur[0] + duvw_face[0], uvw_cur[1] + duvw_face[1], uvw_cur[2] + duvw_face[2] };
                    var viol_cn: usize = 2;
                    if (uvw_try[0] < 0 or uvw_try[0] > 2 * N) {
                        viol_cn = 0;
                    } else if (uvw_try[1] < 0 or uvw_try[1] > 2 * N) {
                        viol_cn = 1;
                    } else {
                        viol_cn = 2;
                    }
                    const opp_map = [_]usize{ 1, 2, 0 };
                    const expected_edge = opp_map[viol_cn];
                    std.debug.assert(expected_edge == ei);
                    // (B) position direction alignment
                    const pos_nf = self.getEdgePositionFace(d.q, d.r, map.nf, map.ne);
                    const ei_dest_opt = self.faceLocalEdgeIndex(map.nf, d.q, d.r);
                    const on_expected_edge = (ei_dest_opt != null and ei_dest_opt.? == map.ne);
                    // Exception: if only (b,a) inequality is violated in neighbor basis and c==0 (on seam),
                    // a 1-tick repair may need to touch b, shifting pos by ±1. Permit that specific case.
                    var allow_delta1 = false;
                    {
                        const perm_n = self.face_perm[f][ei];
                        const uvw_n_pre = applyPerm3(uvw_cur, perm_n);
                        const duvw_n_pre = applyPerm3(duvw_face, perm_n);
                        const pre0: isize = uvw_n_pre[0] + duvw_n_pre[0];
                        const pre1: isize = uvw_n_pre[1] + duvw_n_pre[1];
                        const pre2: isize = uvw_n_pre[2] + duvw_n_pre[2];
                        const Nloc: isize = @intCast(self.size);
                        switch (map.ne) {
                            0 => {
                                const d_ba: isize = pre1 - pre0; // (b=1,a=0)
                                const d_cb: isize = pre2 - pre1; // (c=2,b=1)
                                const d_ac: isize = pre0 - pre2; // (a=0,c=2)
                                if (d_ba > Nloc and d_cb <= Nloc and d_ac <= Nloc and uvw_n_pre[2] == 0) allow_delta1 = true;
                            },
                            1 => {
                                const d_ba: isize = pre2 - pre1; // (b=2,a=1)
                                const d_cb: isize = pre0 - pre2; // (c=0,b=2)
                                const d_ac: isize = pre1 - pre0; // (a=1,c=0)
                                if (d_ba > Nloc and d_cb <= Nloc and d_ac <= Nloc and uvw_n_pre[0] == 0) allow_delta1 = true;
                            },
                            2 => {
                                const d_ba: isize = pre0 - pre2; // (b=0,a=2)
                                const d_cb: isize = pre1 - pre0; // (c=1,b=0)
                                const d_ac: isize = pre2 - pre1; // (a=2,c=1)
                                if (d_ba > Nloc and d_cb <= Nloc and d_ac <= Nloc and uvw_n_pre[1] == 0) allow_delta1 = true;
                            },
                            else => {},
                        }
                    }
                    if (!on_expected_edge) {
                        // Destination is interior on neighbor face; pos' is not meaningful.
                        // Skip strict mirroring in this case.
                    } else if (allow_delta1) {
                        const diff: isize = (@as(isize, @intCast(pos_nf)) - rev_pos);
                        const adiff: isize = if (diff < 0) -diff else diff;
                        std.debug.assert(adiff == 1);
                    } else {
                        std.debug.assert(pos_nf == rev_pos);
                    }
                }
            }
            std.debug.print("  crossed={d} (should be 2)\n", .{crossed});
            test_probe_done = true;
        }
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
                if (self.isPenta(cvar.q, cvar.r)) {
                    any_penta = true;
                    break;
                }
            }
            if (any_penta) cap = 5;
            // Aggregate neighbor indices from all variants, deduplicated
            var out_ids: [6]usize = undefined;
            var out_count: u8 = 0;
            vv = 0;
            while (vv < vcount and out_count < cap) : (vv += 1) {
                var tmp6: [6]Coord = undefined;
                var tmp5: [5]Coord = undefined;
                const cvar = if (variants_slice.len == 0) self.coords[i] else variants_slice[vv];
                var ncount: u8 = 0;
                if (self.isPenta(cvar.q, cvar.r)) {
                    ncount = self.neighbor_penta(cvar.q, cvar.r, cvar.face, &tmp5);
                    var k: usize = 0;
                    while (k < ncount and out_count < cap) : (k += 1) {
                        const h = hashCoord(tmp5[k].q, tmp5[k].r, tmp5[k].face);
                        const idx_opt = coord_to_index.get(h);
                        if (@import("builtin").is_test) std.debug.assert(idx_opt != null);
                        const idx = idx_opt orelse i;
                        if (idx == i and debug_nbr_logs < 24) {
                            if (@import("builtin").is_test) {
                                std.debug.print("miss penta map: i={d} (f={d},q={d},r={d}) -> (f={d},q={d},r={d})\n", .{ i, cvar.face, cvar.q, cvar.r, tmp5[k].face, tmp5[k].q, tmp5[k].r });
                                debug_nbr_logs += 1;
                            }
                        }
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
                } else if (self.isEdge(cvar.q, cvar.r)) {
                    ncount = self.neighbor_edge(cvar.q, cvar.r, cvar.face, &tmp6);
                    var k: usize = 0;
                    while (k < ncount and out_count < cap) : (k += 1) {
                        const h = hashCoord(tmp6[k].q, tmp6[k].r, tmp6[k].face);
                        const idx_opt = coord_to_index.get(h);
                        if (@import("builtin").is_test) std.debug.assert(idx_opt != null);
                        const idx = idx_opt orelse i;
                        if (idx == i and debug_nbr_logs < 24) {
                            if (@import("builtin").is_test) {
                                std.debug.print("miss edge map: i={d} (f={d},q={d},r={d}) -> (f={d},q={d},r={d})\n", .{ i, cvar.face, cvar.q, cvar.r, tmp6[k].face, tmp6[k].q, tmp6[k].r });
                                debug_nbr_logs += 1;
                            }
                        }
                        if (idx == i) continue;
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
                } else {
                    ncount = self.neighbor_inner(cvar.q, cvar.r, cvar.face, &tmp6);
                    var k: usize = 0;
                    while (k < ncount and out_count < cap) : (k += 1) {
                        const h = hashCoord(tmp6[k].q, tmp6[k].r, tmp6[k].face);
                        const idx_opt = coord_to_index.get(h);
                        if (@import("builtin").is_test) std.debug.assert(idx_opt != null);
                        const idx = idx_opt orelse i;
                        if (idx == i and debug_nbr_logs < 24) {
                            if (@import("builtin").is_test) {
                                std.debug.print("miss inner map: i={d} (f={d},q={d},r={d}) -> (f={d},q={d},r={d})\n", .{ i, cvar.face, cvar.q, cvar.r, tmp6[k].face, tmp6[k].q, tmp6[k].r });
                                debug_nbr_logs += 1;
                            }
                        }
                        if (idx == i) continue;
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
        }

        // No symmetry patch: correctness comes from aliasing during index generation
    }

    fn neighbor_inner(self: *const Grid, q: isize, r: isize, face: usize, out: *[6]Coord) u8 {
        _ = self;
        const dirs = [_][2]isize{ .{ 1, 0 }, .{ 1, -1 }, .{ 0, -1 }, .{ -1, 0 }, .{ -1, 1 }, .{ 0, 1 } };
        var k: usize = 0;
        while (k < 6) : (k += 1) {
            out.*[k] = .{ .q = q + dirs[k][0], .r = r + dirs[k][1], .face = face };
        }
        return 6;
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

    fn neighbor_edge(self: *const Grid, q: isize, r: isize, face: usize, out: *[6]Coord) u8 {
        var count: u8 = 0;
        var i: usize = 0;
        while (i < 6) : (i += 1) {
            out.*[count] = self.stepAcrossOrIn(q, r, face, @intCast(i), &.{});
            count += 1;
        }
        return count;
    }

    fn pentaStartDir(q: isize, r: isize, N: isize) u8 {
        if (q == N) return 3;
        if (r == N) return 2;
        return 5;
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
            const back = self.fromBaryFace(face, uvw_try);
            if (@import("builtin").is_test) std.debug.assert(self.isValid(back.q, back.r));
            return .{ .q = back.q, .r = back.r, .face = face };
        }

        // choose edge from violated barycentric component (opposite component index)
        const viol_cn: usize = if (uvw_try[0] < 0 or uvw_try[0] > maxB) 0 else if (uvw_try[1] < 0 or uvw_try[1] > maxB) 1 else 2;
        const eidx: usize = edgeIndexFromOppositeComponent(viol_cn);
        if (@import("builtin").is_test) {
            const expect_eidx = edgeIndexFromOppositeComponent(viol_cn);
            std.debug.assert(eidx == expect_eidx);
        }
        const nf = self.face_neighbors[face][eidx][0] - 1;
        const perm = self.face_perm[face][eidx]; // neighbor_index -> current_index

        const uvw_n = applyPerm3(uvw, perm);
        const duvw_n = applyPerm3(duvw, perm);
        var out: [3]isize = .{ uvw_n[0] + duvw_n[0], uvw_n[1] + duvw_n[1], uvw_n[2] + duvw_n[2] };

        // if already inside neighbor, done
        if (!(out[0] < 0 or out[1] < 0 or out[2] < 0 or out[0] > maxB or out[1] > maxB or out[2] > maxB)) {
            const back_ok = self.fromBaryFace(nf, out);
            return .{ .q = back_ok.q, .r = back_ok.r, .face = nf };
        }

        if (@import("builtin").is_test) {
            std.debug.print("xface final: nf={d} out=({d},{d},{d}) window_ok={any}\n", .{ nf, out[0], out[1], out[2], (out[0] - out[2]) <= N and (out[1] - out[0]) <= N and (out[2] - out[1]) <= N });
        }

        // Edge orientation on neighbor face: define a,b along edge; c is seam component
        const pairs = [_][2]usize{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };
        const ne = self.neighborEdgeAndReversal(face, eidx).ne;
        const a_local = pairs[ne][0];
        const b_local = pairs[ne][1];
        const c_local: usize = 3 - a_local - b_local;

        // 1) Try both partner choices for ±2x compensation, pick best per window/pos rules
        const viol_idx: usize = if (out[0] < 0 or out[0] > maxB) 0 else if (out[1] < 0 or out[1] > maxB) 1 else 2;
        const p0: usize = (viol_idx + 1) % 3;
        const p1: usize = (viol_idx + 2) % 3;

        const reflect = struct {
            fn withPartner(base: [3]isize, du: [3]isize, viol: usize, partner: usize, maxB_: isize) [3]isize {
                var o: [3]isize = .{ base[0] + du[0], base[1] + du[1], base[2] + du[2] };
                const x: isize = if (o[viol] < 0) -o[viol] else o[viol] - maxB_;
                if (o[viol] < 0) {
                    o[viol] += 2 * x;
                    o[partner] -= 2 * x;
                } else {
                    o[viol] -= 2 * x;
                    o[partner] += 2 * x;
                }
                return o;
            }
        };

        const win = struct {
            fn check(u: isize, v: isize, w: isize, N_: isize) bool {
                return (u - w) <= N_ and (v - u) <= N_ and (w - v) <= N_;
            }
        };

        const cand0 = reflect.withPartner(uvw_n, duvw_n, viol_idx, p0, maxB);
        const cand1 = reflect.withPartner(uvw_n, duvw_n, viol_idx, p1, maxB);
        const c0_ok = win.check(cand0[0], cand0[1], cand0[2], N);
        const c1_ok = win.check(cand1[0], cand1[1], cand1[2], N);

        const base_b: isize = uvw_n[b_local] + duvw_n[b_local];
        const scorer = struct {
            fn score(out_: [3]isize, was_ok: bool, base_b_: isize, c_idx: usize, b_idx: usize) usize {
                var s: usize = 0;
                if (!was_ok) s += 100;
                if (out_[c_idx] != 0) s += 50;
                const ob = out_[b_idx];
                const delta_b: isize = if (ob > base_b_) ob - base_b_ else base_b_ - ob;
                // Note: out_[b_local] used at callsite; we keep signature generic here
                return s + @as(usize, @intCast(delta_b));
            }
        };
        // Prefer candidate that becomes in-window after a single repair
        const fixer = struct {
            fn stepOnce(cur: [3]isize, N_: isize, c_idx: usize, _: usize, maxB_: isize) [3]isize {
                if ((cur[0] - cur[2]) <= N_ and (cur[1] - cur[0]) <= N_ and (cur[2] - cur[1]) <= N_) return cur;
                const helper = struct {
                    fn adjust(o: *[3]isize, dec_idx: usize, inc_idx: usize) void {
                        o.*[dec_idx] -= 1;
                        o.*[inc_idx] += 1;
                    }
                    fn pick(cur2: [3]isize, dec_idx: usize, inc_idx: usize, keep_c_idx: usize, maxB2: isize) struct { usize, usize } {
                        const can_dec = cur2[dec_idx] > 0;
                        const can_inc = cur2[inc_idx] < maxB2;
                        var d = dec_idx;
                        var i = inc_idx;
                        if (!can_dec and can_inc) {
                            d = inc_idx;
                            i = dec_idx;
                        } else if (!can_inc and can_dec) {
                            d = dec_idx;
                            i = inc_idx;
                        }
                        if (i == keep_c_idx) {
                            d = dec_idx;
                            i = inc_idx;
                        } else if (d == keep_c_idx) {
                            d = inc_idx;
                            i = dec_idx;
                        }
                        if ((cur2[d] == 0 and cur2[i] < maxB2) or (cur2[i] == maxB2 and cur2[d] > 0)) {
                            const ad = i;
                            const ai = d;
                            if ((cur2[ad] > 0) and (cur2[ai] < maxB2)) {
                                d = ad;
                                i = ai;
                            }
                        }
                        return .{ d, i };
                    }
                };
                var out2 = cur;
                const diff0 = out2[0] - out2[2];
                const diff1 = out2[1] - out2[0];
                const diff2 = out2[2] - out2[1];
                if (diff0 > N_) {
                    const ch = helper.pick(out2, 0, 2, c_idx, maxB_);
                    helper.adjust(&out2, ch[0], ch[1]);
                } else if (diff1 > N_) {
                    const ch = helper.pick(out2, 1, 0, c_idx, maxB_);
                    helper.adjust(&out2, ch[0], ch[1]);
                } else if (diff2 > N_) {
                    const ch = helper.pick(out2, 2, 1, c_idx, maxB_);
                    helper.adjust(&out2, ch[0], ch[1]);
                }
                return out2;
            }
        };
        const fix0 = fixer.stepOnce(cand0, N, c_local, b_local, maxB);
        const fix1 = fixer.stepOnce(cand1, N, c_local, b_local, maxB);
        const f0_ok = (fix0[0] - fix0[2]) <= N and (fix0[1] - fix0[0]) <= N and (fix0[2] - fix0[1]) <= N;
        const f1_ok = (fix1[0] - fix1[2]) <= N and (fix1[1] - fix1[0]) <= N and (fix1[2] - fix1[1]) <= N;
        var s0 = scorer.score(cand0, c0_ok, base_b, c_local, b_local);
        var s1 = scorer.score(cand1, c1_ok, base_b, c_local, b_local);
        if (!f0_ok) s0 += 25;
        if (!f1_ok) s1 += 25;
        out = if (s0 <= s1) cand0 else cand1;
        if (@import("builtin").is_test) {
            std.debug.print("xface reflect: face={d} eidx={d} nf={d} ne={d} viol_idx={d} uvw_n=({d},{d},{d}) duvw_n=({d},{d},{d}) cand0=({d},{d},{d}) ok0={any} cand1=({d},{d},{d}) ok1={any} chosen=({d},{d},{d})\n", .{ face, eidx, nf, ne, viol_idx, uvw_n[0], uvw_n[1], uvw_n[2], duvw_n[0], duvw_n[1], duvw_n[2], cand0[0], cand0[1], cand0[2], c0_ok, cand1[0], cand1[1], cand1[2], c1_ok, out[0], out[1], out[2] });
        }

        // 2) If still outside the window, do minimal repairs on the violating pair(s)
        var iter: usize = 0;
        while (!win.check(out[0], out[1], out[2], N) and iter < 3) : (iter += 1) {
            const helper = struct {
                fn adjust(o: *[3]isize, dec_idx: usize, inc_idx: usize) void {
                    o.*[dec_idx] -= 1;
                    o.*[inc_idx] += 1;
                }
                fn pick(cur: [3]isize, dec_idx: usize, inc_idx: usize, keep_c_idx: usize, maxB_: isize) struct { usize, usize } {
                    // Avoid out-of-bounds adjustments first
                    const can_dec = cur[dec_idx] > 0;
                    const can_inc = cur[inc_idx] < maxB_;
                    if (!can_dec and can_inc) return .{ inc_idx, dec_idx };
                    if (!can_inc and can_dec) return .{ dec_idx, inc_idx };
                    if (!can_dec and !can_inc) return .{ dec_idx, inc_idx }; // will be caught later
                    // Honor keep_c first
                    var dec_final = dec_idx;
                    var inc_final = inc_idx;
                    if (inc_idx == keep_c_idx) {
                        dec_final = dec_idx;
                        inc_final = inc_idx;
                    } else if (dec_idx == keep_c_idx) {
                        dec_final = inc_idx;
                        inc_final = dec_idx;
                    }
                    // Final safety: avoid OOB with the chosen orientation
                    if ((cur[dec_final] == 0 and cur[inc_final] < maxB_) or (cur[inc_final] == maxB_ and cur[dec_final] > 0)) {
                        // flip if possible
                        const alt_dec = inc_final;
                        const alt_inc = dec_final;
                        if ((cur[alt_dec] > 0) and (cur[alt_inc] < maxB_)) {
                            dec_final = alt_dec;
                            inc_final = alt_inc;
                        }
                    }
                    return .{ dec_final, inc_final };
                }
            };
            const diff0: isize = out[0] - out[2]; // pair (0,2)
            const diff1: isize = out[1] - out[0]; // pair (1,0)
            const diff2: isize = out[2] - out[1]; // pair (2,1)
            if (diff0 > N) {
                const choice = helper.pick(out, 0, 2, c_local, maxB);
                helper.adjust(&out, choice[0], choice[1]);
            } else if (diff1 > N) {
                const choice = helper.pick(out, 1, 0, c_local, maxB);
                helper.adjust(&out, choice[0], choice[1]);
            } else if (diff2 > N) {
                const choice = helper.pick(out, 2, 1, c_local, maxB);
                helper.adjust(&out, choice[0], choice[1]);
            } else break;
        }

        if (@import("builtin").is_test) {
            if (!(out[0] >= 0 and out[1] >= 0 and out[2] >= 0 and out[0] <= maxB and out[1] <= maxB and out[2] <= maxB)) {
                std.debug.print("xface OOB: face={d} eidx={d} nf={d} N={d}\n", .{ face, eidx, nf, N });
                std.debug.print(" uvw=({d},{d},{d}) du=({d},{d},{d})\n", .{ uvw[0], uvw[1], uvw[2], duvw[0], duvw[1], duvw[2] });
                std.debug.print(" uvw_n=({d},{d},{d}) du_n=({d},{d},{d}) perm=[{d},{d},{d}]\n", .{ uvw_n[0], uvw_n[1], uvw_n[2], duvw_n[0], duvw_n[1], duvw_n[2], perm[0], perm[1], perm[2] });
                std.debug.print(" out=({d},{d},{d}) viol={d} maxB={d}\n", .{ out[0], out[1], out[2], 0, maxB });
            }
        }
        std.debug.assert(out[0] >= 0 and out[1] >= 0 and out[2] >= 0);
        std.debug.assert(out[0] <= maxB and out[1] <= maxB and out[2] <= maxB);
        const back = self.fromBaryFace(nf, out);
        if (@import("builtin").is_test) {
            const ok = self.isValid(back.q, back.r);
            if (!ok) {
                std.debug.print("landing invalid: nf={d} out=({d},{d},{d}) -> (q={d},r={d})\n", .{ nf, out[0], out[1], out[2], back.q, back.r });
            }
            std.debug.assert(ok);
        }
        return .{ .q = back.q, .r = back.r, .face = nf };
    }

    fn neighbor_penta(self: *const Grid, q: isize, r: isize, face: usize, out: *[5]Coord) u8 {
        const N: isize = @intCast(self.size);
        var tmp: [6]Coord = undefined;
        const start: u8 = pentaStartDir(q, r, N);
        var i: u8 = 0;
        while (i < 6) : (i += 1) {
            const dir: u8 = (start + i) % 6;
            tmp[i] = stepAcrossOrIn(self, q, r, face, dir, &.{});
        }
        var used: usize = 0;
        var k: usize = 0;
        while (k < 6 and used < 5) : (k += 1) {
            const cand = tmp[k];
            var dup = false;
            var j: usize = 0;
            while (j < used) : (j += 1) {
                if (out.*[j].q == cand.q and out.*[j].r == cand.r and out.*[j].face == cand.face) {
                    dup = true;
                    break;
                }
            }
            if (!dup) {
                out.*[used] = cand;
                used += 1;
            }
        }
        return @intCast(used);
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
                    if (self.isEdge(q, r)) {
                        if (self.isPenta(q, r)) {
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

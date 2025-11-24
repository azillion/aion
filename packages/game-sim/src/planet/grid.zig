const std = @import("std");
const Vec3 = struct { x: f64, y: f64, z: f64 };
const QR = struct { q: isize, r: isize };
pub const Coord = struct { q: isize, r: isize, face: usize };
pub const NeighborSet = struct { ids: [6]usize, count: u8 };

inline fn legacyRemoved(comptime msg: []const u8) noreturn {
    @panic(msg);
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

pub const Grid = struct {
    allocator: std.mem.Allocator,
    size: usize,
    tile_count: usize,
    coords: []Coord,
    neighbors: []NeighborSet,
    vertices: ?[]f32,
    elevations: ?[]f32,
    indices: ?[]u32,

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
        return Grid{
            .allocator = allocator,
            .size = size,
            .tile_count = 0,
            .coords = &[_]Coord{},
            .neighbors = &[_]NeighborSet{},
            .vertices = null,
            .elevations = null,
            .indices = null,
        };
    }

    pub fn deinit(self: *Grid) void {
        if (self.coords.len > 0) self.allocator.free(self.coords);
        if (self.neighbors.len > 0) self.allocator.free(self.neighbors);
        if (self.vertices) |v| self.allocator.free(v);
        if (self.elevations) |e| self.allocator.free(e);
        if (self.indices) |i| self.allocator.free(i);
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

    inline fn toBary(N: isize, q: isize, r: isize) [3]isize {
        const u = q + N;
        const v = r + N;
        const w = N - q - r;
        return .{ u, v, w };
    }

    pub inline fn toBary_public(N: isize, q: isize, r: isize) [3]isize {
        return toBary(N, q, r);
    }

    inline fn fromBary(N: isize, uvw: [3]isize) QR {
        return QR{ .q = uvw[0] - N, .r = uvw[1] - N };
    }

    pub fn faceCenterPosition(self: *const Grid, q: isize, r: isize, face: usize) Vec3 {
        const N: isize = @intCast(self.size);
        const uvw = toBary(N, q, r);
        const tri = icosa_face_vertices[face];
        const a = icosa_vertices[tri[0] - 1];
        const b = icosa_vertices[tri[1] - 1];
        const c = icosa_vertices[tri[2] - 1];
        const denom = @as(f64, @floatFromInt(3 * N));
        const wa = @as(f64, @floatFromInt(uvw[0])) / denom;
        const wb = @as(f64, @floatFromInt(uvw[1])) / denom;
        const wc = @as(f64, @floatFromInt(uvw[2])) / denom;
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

    pub fn getGlobalVertex(self: *const Grid, face: usize, q: isize, r: isize) ?usize {
        const p = self.faceCenterPosition(q, r, face);
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

    pub fn isValid(self: *const Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const N = @as(isize, @intCast(self.size));
        return iabs(q) <= N and iabs(r) <= N and iabs(s) <= N;
    }

    pub fn isEdge(self: *const Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const size_i = @as(isize, @intCast(self.size));
        return q - s == size_i or r - q == size_i or s - r == size_i;
    }

    pub fn isPentaAt(self: *const Grid, face: usize, q: isize, r: isize) bool {
        return self.isPentaFace(face, q, r);
    }

    fn isPentaFace(self: *const Grid, _: usize, q: isize, r: isize) bool {
        const N: isize = @intCast(self.size);
        const uvw = toBary(N, q, r);
        var zero_count: usize = 0;
        var has_2N = false;
        var has_N = false;
        var i: usize = 0;
        while (i < 3) : (i += 1) {
            if (uvw[i] == 0) zero_count += 1 else if (uvw[i] == 2 * N) has_2N = true else if (uvw[i] == N) has_N = true;
        }
        return zero_count == 1 and has_2N and has_N;
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

    pub fn hashCoord(q: isize, r: isize, face: usize) u64 {
        const sz: i64 = 1 << 20;
        const max_coord: isize = @intCast(sz);
        std.debug.assert(q >= -max_coord and q < max_coord);
        std.debug.assert(r >= -max_coord and r < max_coord);
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

    inline fn iabs(x: isize) isize {
        return if (x < 0) -x else x;
    }
};

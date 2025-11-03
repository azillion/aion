const std = @import("std");
const Vec3 = struct { x: f64, y: f64, z: f64 };
const UNSET = std.math.maxInt(usize);
const TransformResult = struct { neighbor_face: usize, neighbor_weights: [3]f64 };
pub const Coord = struct { q: isize, r: isize, face: usize };
pub const NeighborSet = struct { ids: [6]usize, count: u8 };
comptime {
    // 22 (face) + 21 (q) + 21 (r) = 64
    std.debug.assert((22 + 21 + 21) == 64);
}

pub const icosa_vertices = [_]Vec3{
    .{ .x = -0.607026219367981, .y = 0.794681787490845, .z = 0.00008639693260193 },
    .{ .x = 0.303695321083069, .y = 0.794681787490845, .z = 0.525594890117645 },
    .{ .x = -0.303695321083069, .y = -0.794681787490845, .z = -0.525594890117645 },
    .{ .x = 0.607026219367981, .y = -0.794681787490845, .z = -0.00008639693260193 },
    .{ .x = 0.303438216447830, .y = 0.794599771499634, .z = -0.525867342948914 },
    .{ .x = 0.490907043218613, .y = -0.187680929899216, .z = -0.850756168365479 },
    .{ .x = -0.490907043218613, .y = 0.187680929899216, .z = 0.850756168365479 },
    .{ .x = -0.303438216447830, .y = -0.794599771499634, .z = 0.525867342948914 },
    .{ .x = -0.491322875022888, .y = 0.187548249959946, .z = -0.850545346736908 },
    .{ .x = 0.982255339622498, .y = 0.187548249959946, .z = -0.00025475025177002 },
    .{ .x = -0.982255339622498, .y = -0.187548249959946, .z = 0.00025475025177002 },
    .{ .x = 0.491322875022888, .y = -0.187548249959946, .z = 0.850545346736908 },
};

pub const icosa_face_vertices = [_][3]usize{ .{ 1, 2, 5 }, .{ 2, 10, 5 }, .{ 5, 10, 6 }, .{ 6, 10, 4 }, .{ 3, 4, 8 }, .{ 4, 3, 6 }, .{ 8, 11, 3 }, .{ 1, 9, 11 }, .{ 1, 5, 9 }, .{ 9, 3, 11 }, .{ 9, 5, 6 }, .{ 9, 6, 3 }, .{ 2, 1, 7 }, .{ 12, 2, 7 }, .{ 4, 10, 12 }, .{ 7, 11, 8 }, .{ 4, 12, 8 }, .{ 12, 7, 8 }, .{ 7, 1, 11 }, .{ 10, 2, 12 } };

pub const icosa_face_neighbors = [_][3][2]usize{ .{ .{ 2, 3 }, .{ 13, 2 }, .{ 9, 2 } }, .{ .{ 3, 2 }, .{ 20, 2 }, .{ 1, 1 } }, .{ .{ 4, 3 }, .{ 1, 3 }, .{ 11, 3 } }, .{ .{ 5, 3 }, .{ 2, 3 }, .{ 15, 3 } }, .{ .{ 6, 3 }, .{ 7, 2 }, .{ 12, 3 } }, .{ .{ 7, 1 }, .{ 4, 1 }, .{ 11, 1 } }, .{ .{ 8, 3 }, .{ 16, 3 }, .{ 6, 3 } }, .{ .{ 9, 3 }, .{ 19, 3 }, .{ 1, 2 } }, .{ .{ 10, 3 }, .{ 1, 1 }, .{ 11, 2 } }, .{ .{ 11, 1 }, .{ 8, 1 }, .{ 12, 2 } }, .{ .{ 12, 1 }, .{ 3, 1 }, .{ 10, 1 } }, .{ .{ 6, 1 }, .{ 5, 1 }, .{ 10, 2 } }, .{ .{ 14, 3 }, .{ 1, 3 }, .{ 19, 1 } }, .{ .{ 15, 3 }, .{ 2, 1 }, .{ 13, 1 } }, .{ .{ 16, 3 }, .{ 4, 2 }, .{ 20, 3 } }, .{ .{ 17, 3 }, .{ 19, 2 }, .{ 8, 2 } }, .{ .{ 18, 3 }, .{ 15, 2 }, .{ 4, 1 } }, .{ .{ 14, 1 }, .{ 17, 1 }, .{ 12, 1 } }, .{ .{ 13, 2 }, .{ 8, 1 }, .{ 17, 2 } }, .{ .{ 16, 2 }, .{ 14, 2 }, .{ 2, 2 } } };

// Canonical axial direction list (CW): (+q), (+r), (-q,+s), (-q), (-r), (+q,-s)
const DIRS = [_][2]isize{ .{ 1, 0 }, .{ 1, -1 }, .{ 0, -1 }, .{ -1, 0 }, .{ -1, 1 }, .{ 0, 1 } };

// Maps a direction index to its reflected index across each of the 3 edges
const DIR_REFLECTION = [_][6]u8{
    .{ 0, 5, 4, 3, 2, 1 }, // Edge 0 (q - s)
    .{ 2, 1, 0, 5, 4, 3 }, // Edge 1 (r - q)
    .{ 4, 3, 2, 1, 0, 5 }, // Edge 2 (s - r)
};

var debug_edge_logs: usize = 0;
var debug_nbr_logs: usize = 0;

pub const Grid = struct {
    allocator: std.mem.Allocator,
    size: usize,
    tile_count: usize,

    // Final topology tables
    coords: []Coord,
    neighbors: []NeighborSet,

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
        var grid = Grid{
            .allocator = allocator,
            .size = size,
            .tile_count = 0,
            .coords = &[_]Coord{},
            .neighbors = &[_]NeighborSet{},
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

        return grid;
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

    fn faceCenterPosition(self: *Grid, q: isize, r: isize, face: usize) Vec3 {
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
        const ei = self.getEdgeIndex(q, r) orelse return null;
        const edge_len: usize = if (self.size > 0) self.size - 1 else 0;
        if (edge_len == 0) return null;
        const pos = self.getEdgePosition(q, r, ei);
        if (pos >= edge_len) return null;
        const idx_a = self.edge_indices[face][ei][pos];
        const map = neighborEdgeAndReversal(face, ei);
        const rev_pos = if (map.rev) (edge_len - 1) - pos else pos;
        if (map.nf >= 20 or map.ne >= 3 or rev_pos >= self.edge_indices[map.nf][map.ne].len) return null;
        const idx_b = self.edge_indices[map.nf][map.ne][rev_pos];
        return .{ .a = idx_a, .b = idx_b };
    }

    // Test helper: expose neighbor edge mapping for tests
    pub fn testNeighborEdgeMap(face: usize, ei: usize) EdgeMap {
        const map = neighborEdgeAndReversal(face, ei);
        return .{ .nf = map.nf, .ne = map.ne, .rev = map.rev };
    }

    // Test helper: expose static neighbor lookup
    pub fn faceNeighborFace(face: usize, ei: usize) usize {
        return icosa_face_neighbors[face][ei][0] - 1;
    }

    // Position 0..(N-2) along interior of an edge (exclude pentagon corners)
    fn getEdgePosition(self: *Grid, q: isize, r: isize, edge: usize) usize {
        std.debug.assert(self.isEdge(q, r));
        std.debug.assert(!self.isPenta(q, r));
        const N: isize = @intCast(self.size);
        const pos: isize = switch (edge) {
            // Edge 0: q - s == N -> q in [1..N-1]
            0 => q - 1,
            // Edge 1: r - q == N -> q in [-(N-1)..-1]
            1 => (-q) - 1,
            // Edge 2: s - r == N -> r in [-(N-1)..-1]
            2 => (-r) - 1,
            else => 0,
        };
        std.debug.assert(pos >= 0 and pos <= (N - 2));
        return @intCast(pos);
    }

    // getVertexIndex no longer used; kept for compatibility if referenced elsewhere
    fn getVertexIndex(self: *Grid, q: isize, r: isize) usize {
        const N: isize = @intCast(self.size);
        if (q == N) return 0;
        if (r == N) return 1;
        return 2;
    }

    fn neighborEdgeAndReversal(face: usize, ei: usize) EdgeMap {
        const cur = icosa_face_vertices[face];
        const a_local = (ei + 1) % 3;
        const b_local = (ei + 2) % 3;
        const ga = cur[a_local];
        const gb = cur[b_local];

        const nf = icosa_face_neighbors[face][ei][0] - 1;
        const nfv = icosa_face_vertices[nf];

        var ia: usize = 3;
        var ib: usize = 3;
        var i: usize = 0;
        while (i < 3) : (i += 1) {
            if (nfv[i] == ga) ia = i;
            if (nfv[i] == gb) ib = i;
        }
        // Fallback if shared vertices not found (shouldn't happen)
        if (ia == 3 or ib == 3) {
            const ne_tbl = icosa_face_neighbors[face][ei][1] - 1;
            return .{ .nf = nf, .ne = ne_tbl, .rev = true };
        }
        var ne: usize = 0;
        while (ne < 3) : (ne += 1) {
            if (ne != ia and ne != ib) break;
        }
        const next = [_]usize{ 1, 2, 0 };
        const ia_next = next[ia];
        const rev = !(ia_next == ib);
        return .{ .nf = nf, .ne = ne, .rev = rev };
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
        const ei = self.getEdgeIndex(q, r) orelse return .{ current_index, current_index + 1 };
        const edge_len: usize = if (self.size > 0) self.size - 1 else 0;
        if (edge_len == 0) return .{ current_index, current_index + 1 };
        const pos = self.getEdgePosition(q, r, ei);
        if (pos >= edge_len) return .{ current_index, current_index + 1 };
        if (self.edge_indices[face][ei][pos] != UNSET) return .{ self.edge_indices[face][ei][pos], current_index };
        self.edge_indices[face][ei][pos] = current_index;
        const map = neighborEdgeAndReversal(face, ei);
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
        var i: usize = 0;
        while (i < self.tile_count) : (i += 1) {
            // initialize ring to self and reset count to 0 before filling
            @memset(&self.neighbors[i].ids, i);
            self.neighbors[i].count = 0;
            // Collect all coordinate variants that alias to this index (handles edge/penta dedup)
            var variants: [8]Coord = undefined;
            var vcount: usize = 0;
            var itv = coord_to_index.iterator();
            while (itv.next()) |e| {
                if (e.value_ptr.* == i) {
                    const cvar = unhashCoord(e.key_ptr.*);
                    if (vcount < variants.len) {
                        variants[vcount] = .{ .q = cvar.q, .r = cvar.r, .face = cvar.face };
                        vcount += 1;
                    }
                }
            }
            if (vcount == 0) {
                variants[0] = self.coords[i];
                vcount = 1;
            }
            // Determine capacity: if any variant is a penta, cap=5 else 6
            var cap: u8 = 6;
            var any_penta = false;
            var vv: usize = 0;
            while (vv < vcount) : (vv += 1) {
                if (self.isPenta(variants[vv].q, variants[vv].r)) {
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
                const cvar = variants[vv];
                var ncount: u8 = 0;
                if (self.isPenta(cvar.q, cvar.r)) {
                    ncount = self.neighbor_penta(cvar.q, cvar.r, cvar.face, &tmp5);
                    var k: usize = 0;
                    while (k < ncount and out_count < cap) : (k += 1) {
                        const h = hashCoord(tmp5[k].q, tmp5[k].r, tmp5[k].face);
                        const idx = coord_to_index.get(h) orelse i;
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
                        const idx = coord_to_index.get(h) orelse i;
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
                        const idx = coord_to_index.get(h) orelse i;
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
            // Write results and pad
            var w: usize = 0;
            while (w < out_count) : (w += 1) self.neighbors[i].ids[w] = out_ids[w];
            while (w < cap) : (w += 1) self.neighbors[i].ids[w] = i;
            self.neighbors[i].count = cap;
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
        const nq = q + dq;
        const nr = r + dr;
        if (self.isValid(nq, nr)) return .{ .q = nq, .r = nr, .face = face };

        var eidx: usize = 0;
        {
            const N: isize = @intCast(self.size);
            const s = -(q + r);
            var minv = q + N + dq;
            if (r + N + dr < minv) {
                minv = r + N + dr;
                eidx = 1;
            }
            if (s + N - (dq + dr) < minv) {
                eidx = 2;
            }
        }

        const mapped_origin = transformCoord(self, q, r, eidx);
        const Step = struct { q: isize, r: isize };
        const mapped_step: Step = switch (eidx) {
            0 => Step{ .q = dq + dr, .r = -dr },
            1 => Step{ .q = -dq, .r = dq + dr },
            2 => Step{ .q = -dr, .r = -dq },
            else => unreachable,
        };
        const nf = icosa_face_neighbors[face][eidx][0] - 1;
        return .{ .q = mapped_origin.q + mapped_step.q, .r = mapped_origin.r + mapped_step.r, .face = nf };
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
            if (ring.count < 2) continue;
            var j: u8 = 0;
            while (j < ring.count) : (j += 1) {
                const b_idx = ring.ids[j];
                const c_idx = ring.ids[@intCast((@as(u16, j) + 1) % @as(u16, ring.count))];
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

const std = @import("std");
const math = @import("../math");
const Vec3 = math.Vec3;

pub const icosa_vertices = [_]Vec3{
    Vec3{ .x = -0.60702621936798100, .y = 0.79468178749084500, .z = 0.00008639693260193 },
    Vec3{ .x = 0.30369532108306900, .y = 0.79468178749084500, .z = 0.52559489011764500 },
    Vec3{ .x = -0.30369532108306900, .y = 0.79468178749084500, .z = -0.52559489011764500 },
    Vec3{ .x = 0.60702621936798100, .y = 0.79468178749084500, .z = -0.00008639693260193 },
    Vec3{ .x = -0.30369532108306900, .y = 0.79468178749084500, .z = 0.52559489011764500 },
    Vec3{ .x = 0.00008639693260193, .y = 0.79468178749084500, .z = 0.60702621936798100 },
    Vec3{ .x = 0.30369532108306900, .y = 0.79468178749084500, .z = -0.52559489011764500 },
    Vec3{ .x = -0.00008639693260193, .y = 0.79468178749084500, .z = -0.60702621936798100 },
    Vec3{ .x = -0.52559489011764500, .y = 0.30369532108306900, .z = 0.79468178749084500 },
    Vec3{ .x = 0.52559489011764500, .y = 0.30369532108306900, .z = 0.79468178749084500 },
    Vec3{ .x = 0.60702621936798100, .y = 0.00008639693260193, .z = 0.79468178749084500 },
    Vec3{ .x = 0.52559489011764500, .y = -0.30369532108306900, .z = 0.79468178749084500 },
};

pub const icosa_face_vertices = [_][3]usize{ .{ 1, 2, 5 }, .{ 2, 10, 5 }, .{ 5, 10, 6 }, .{ 6, 10, 4 }, .{ 3, 4, 8 }, .{ 4, 3, 6 }, .{ 8, 11, 3 }, .{ 1, 9, 11 }, .{ 1, 5, 9 }, .{ 9, 3, 11 }, .{ 9, 5, 6 }, .{ 9, 6, 3 }, .{ 2, 1, 7 }, .{ 12, 2, 7 }, .{ 4, 10, 12 }, .{ 7, 11, 8 }, .{ 4, 12, 8 }, .{ 12, 7, 8 }, .{ 7, 1, 11 }, .{ 10, 2, 12 } };

pub const icosa_face_neighbors = [_][3][2]usize{ .{ .{ 2, 3 }, .{ 13, 2 }, .{ 9, 2 } }, .{ .{ 3, 2 }, .{ 20, 2 }, .{ 1, 1 } }, .{ .{ 4, 3 }, .{ 1, 3 }, .{ 11, 3 } }, .{ .{ 5, 3 }, .{ 2, 3 }, .{ 15, 3 } }, .{ .{ 6, 3 }, .{ 7, 2 }, .{ 12, 3 } }, .{ .{ 7, 1 }, .{ 4, 1 }, .{ 11, 1 } }, .{ .{ 8, 3 }, .{ 16, 3 }, .{ 6, 3 } }, .{ .{ 9, 3 }, .{ 19, 3 }, .{ 1, 2 } }, .{ .{ 10, 3 }, .{ 1, 1 }, .{ 11, 2 } }, .{ .{ 11, 1 }, .{ 8, 1 }, .{ 12, 2 } }, .{ .{ 12, 1 }, .{ 3, 1 }, .{ 10, 1 } }, .{ .{ 6, 1 }, .{ 5, 1 }, .{ 10, 2 } }, .{ .{ 14, 3 }, .{ 1, 3 }, .{ 19, 1 } }, .{ .{ 15, 3 }, .{ 2, 1 }, .{ 13, 1 } }, .{ .{ 16, 3 }, .{ 4, 2 }, .{ 20, 3 } }, .{ .{ 17, 3 }, .{ 19, 2 }, .{ 8, 2 } }, .{ .{ 18, 3 }, .{ 15, 2 }, .{ 4, 1 } }, .{ .{ 14, 1 }, .{ 17, 1 }, .{ 12, 1 } }, .{ .{ 13, 2 }, .{ 8, 1 }, .{ 17, 2 } }, .{ .{ 16, 2 }, .{ 14, 2 }, .{ 2, 2 } } };

pub const Grid = struct {
    allocator: std.mem.Allocator,
    size: usize,
    tile_count: usize,

    // Data structures for tile information
    neighbors: ?[][]usize,
    coord_map: ?std.hash_map.AutoHashMap(u64, usize),

    // Internal data structures for tracking shared indices during generation
    penta_indices: [12]usize,
    edge_indices: [20][3][]usize,

    pub fn init(allocator: std.mem.Allocator, size: usize) !Grid {
        var grid = Grid{
            .allocator = allocator,
            .size = size,
            .tile_count = 0,
            .neighbors = null,
            .coord_map = null,
            .penta_indices = std.mem.zeroes([12]usize),
            .edge_indices = std.mem.zeroes([20][3][]usize),
        };

        // Initialize data structures
        grid.coord_map = std.hash_map.AutoHashMap(u64, usize).init(allocator);

        // Allocate edge indices arrays
        for (0..20) |face| {
            for (0..3) |edge| {
                const edge_len = size * 2 + 1;
                grid.edge_indices[face][edge] = try allocator.alloc(usize, edge_len);
                @memset(grid.edge_indices[face][edge], std.math.maxInt(usize));
            }
        }

        return grid;
    }

    pub fn deinit(self: *Grid) void {
        if (self.neighbors) |neighbors| {
            for (neighbors) |neighbor_row| {
                self.allocator.free(neighbor_row);
            }
            self.allocator.free(neighbors);
        }
        if (self.coord_map) |*coord_map| {
            coord_map.deinit();
        }
        for (0..20) |face| {
            for (0..3) |edge| {
                self.allocator.free(self.edge_indices[face][edge]);
            }
        }
    }

    // Helper function to hash coordinates for map lookup
    pub fn hashCoord(q: isize, r: isize, face: usize) u64 {
        // Convert signed coordinates to unsigned by adding offset to make them positive
        const offset = 1 << 20; // Large enough offset for our coordinate range
        const q_u = @as(u64, @intCast(q + offset));
        const r_u = @as(u64, @intCast(r + offset));

        return (@as(u64, @intCast(face)) << 44) |
            ((q_u & 0xFFFFF) << 22) |
            (r_u & 0x3FFFFF);
    }

    // Helper function to check if a hex coordinate is valid for a given size
    pub fn isValid(self: *Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const size_i = @as(isize, @intCast(self.size));
        return q - s <= size_i and r - q <= size_i and s - r <= size_i;
    }

    // Helper function to check if a tile is on the main icosahedron edge
    pub fn isEdge(self: *Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const size_i = @as(isize, @intCast(self.size));
        return q - s == size_i or r - q == size_i or s - r == size_i;
    }

    // Helper function to check if a tile is a pentagon corner
    pub fn isPenta(self: *Grid, q: isize, r: isize) bool {
        const s = -(q + r);
        const size_i = @as(isize, @intCast(self.size));
        return q == size_i or r == size_i or s == size_i;
    }

    // Helper function to get which edge a coordinate is on
    fn getEdgeIndex(self: *Grid, q: isize, r: isize) ?usize {
        const size_i = @as(isize, @intCast(self.size));
        const s = -(q + r);

        if (q - s == size_i) return 0; // Edge 0
        if (r - q == size_i) return 1; // Edge 1
        if (s - r == size_i) return 2; // Edge 2

        return null;
    }

    // Helper function to get position along edge
    fn getEdgePosition(self: *Grid, q: isize, r: isize, edge: usize) usize {
        _ = self;
        _ = q;
        _ = r;

        switch (edge) {
            0 => return 0, // Will be calculated properly in full implementation
            1 => return 0,
            2 => return 0,
            else => return 0,
        }
    }

    // Helper function to get vertex index for pentagon
    fn getVertexIndex(self: *Grid, q: isize, r: isize) usize {
        _ = q;
        _ = r;
        _ = self;
        return 0; // Will be implemented properly
    }

    // Resolve index for pentagon tiles (12 special vertices of icosahedron)
    fn resolveIndexForPenta(self: *Grid, q: isize, r: isize, face: usize, current_index: usize) struct { usize, usize } {
        _ = face;

        const vertex_index = self.getVertexIndex(q, r);

        if (self.penta_indices[vertex_index] != std.math.maxInt(usize)) {
            return .{ self.penta_indices[vertex_index], current_index };
        }

        self.penta_indices[vertex_index] = current_index;
        return .{ current_index, current_index + 1 };
    }

    // Resolve index for edge tiles (tiles on the edges between faces)
    fn resolveIndexForEdge(self: *Grid, q: isize, r: isize, face: usize, current_index: usize) struct { usize, usize } {
        const edge_index = self.getEdgeIndex(q, r) orelse return .{ current_index, current_index + 1 };
        const edge_pos = self.getEdgePosition(q, r, edge_index);

        if (self.edge_indices[face][edge_index][edge_pos] != std.math.maxInt(usize)) {
            return .{ self.edge_indices[face][edge_index][edge_pos], current_index };
        }

        self.edge_indices[face][edge_index][edge_pos] = current_index;

        // Also assign to neighboring face's edge
        const neighbor_info = icosa_face_neighbors[face][edge_index];
        const neighbor_face = neighbor_info[0] - 1; // Convert from 1-based to 0-based
        const neighbor_edge = neighbor_info[1] - 1; // Convert from 1-based to 0-based

        if (neighbor_face < 20 and neighbor_edge < 3) {
            const neighbor_pos = self.size * 2 - edge_pos; // Reverse position on neighbor
            self.edge_indices[neighbor_face][neighbor_edge][neighbor_pos] = current_index;
        }

        return .{ current_index, current_index + 1 };
    }

    // Populate indices for all tiles in the grid
    pub fn populateIndices(self: *Grid) !void {
        if (self.coord_map) |*coord_map| {
            var index: usize = 0;
            const size_i = @as(isize, @intCast(self.size));

            // Main loop for assigning indices - iterate from -size to size
            var q = -size_i;
            while (q <= size_i) : (q += 1) {
                var r = -size_i;
                while (r <= size_i) : (r += 1) {
                    // Skip invalid hex coordinates
                    if (!self.isValid(q, r)) continue;

                    // Iterate through all 20 faces
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
                            resolved_index = index;
                            next_index = index + 1;
                        }

                        try coord_map.put(coord_hash, resolved_index);
                        index = next_index;
                    }
                }
            }

            self.tile_count = index;
        }
    }

    // Helper function to get neighbor in a specific direction
    fn getNeighborInDirection(self: *Grid, q: isize, r: isize, face: usize, direction: usize) ?usize {
        _ = self;
        _ = face;

        // Hex grid directions (6 neighbors)
        const directions = [_][2]isize{
            .{ 1, 0 },  .{ 1, -1 }, .{ 0, -1 },
            .{ -1, 0 }, .{ -1, 1 }, .{ 0, 1 },
        };

        if (direction >= directions.len) return null;

        const nq = q + directions[direction][0];
        const nr = r + directions[direction][1];

        // In a full implementation, this would need to handle:
        // 1. Face wrapping
        // 2. Edge cases
        // 3. Coordinate system conversion between faces

        // For now, return a simple hash-based index
        return @as(usize, @intCast(@abs(nq) + @abs(nr) * 1000));
    }

    // Populate neighbors for all tiles
    pub fn populateNeighbors(self: *Grid) !void {
        // Allocate neighbors array
        self.neighbors = try self.allocator.alloc([]usize, self.tile_count);
        for (self.neighbors.?) |*neighbor_row| {
            neighbor_row.* = try self.allocator.alloc(usize, 6); // 6 neighbors per hex
        }

        if (self.coord_map) |*coord_map| {
            var tile_iter = coord_map.iterator();
            while (tile_iter.next()) |entry| {
                const tile_index = entry.key_ptr.*;
                const coord = entry.value_ptr.*;
                const q = coord[0];
                const r = coord[1];
                const face = coord[2];

                // Find all 6 neighbors
                for (0..6) |direction| {
                    if (self.getNeighborInDirection(q, r, @as(usize, @intCast(face)), direction)) |neighbor_index| {
                        self.neighbors.?[tile_index][direction] = neighbor_index;
                    } else {
                        self.neighbors.?[tile_index][direction] = tile_index; // Self as fallback
                    }
                }
            }
        }
    }
};

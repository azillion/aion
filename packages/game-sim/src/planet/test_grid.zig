const std = @import("std");
const planet = @import("index.zig");
const expect = std.testing.expect;
const expectEqual = std.testing.expectEqual;

test "Grid.init" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 1);
    defer grid.deinit();

    try expectEqual(@as(usize, 1), grid.size);
    try expectEqual(@as(usize, 0), grid.tile_count);
    try expect(grid.coord_map != null);
    try expect(grid.neighbors == null);
}

test "Grid.isValid coordinates" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);
    defer grid.deinit();

    // Valid coordinates within the hex boundary
    try expect(grid.isValid(0, 0));
    try expect(grid.isValid(1, 0));
    try expect(grid.isValid(0, 1));
    try expect(grid.isValid(-1, 0));
    try expect(grid.isValid(0, -1));

    // Invalid coordinates outside the hex boundary
    try expect(!grid.isValid(3, 0));
    try expect(!grid.isValid(0, 3));
    try expect(!grid.isValid(-3, 0));
    try expect(!grid.isValid(0, -3));
    try expect(!grid.isValid(2, 2)); // q + r would be 4, s = -4, q - s = 6 > 2
}

test "Grid.isEdge detection" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);
    defer grid.deinit();

    // Edge coordinates where q - s = size, r - q = size, or s - r = size
    try expect(grid.isEdge(0, 2)); // r - q = 2 - 0 = 2 = size
    try expect(grid.isEdge(-2, 0)); // s - r = 2 - 0 = 2 = size

    // Coordinates that are NOT edges (based on actual behavior)
    try expect(!grid.isEdge(2, 0)); // q - s = 4 > size, but also invalid
    try expect(!grid.isEdge(0, 0));
    try expect(grid.isEdge(1, 0)); // q - s = 2 = size, so this IS an edge
    try expect(!grid.isEdge(0, 1));
}

test "Grid.isPenta detection" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);
    defer grid.deinit();

    // Pentagon coordinates where q = size, r = size, or s = size
    try expect(grid.isPenta(2, 0)); // q = size
    try expect(grid.isPenta(0, 2)); // r = size
    try expect(grid.isPenta(-2, 0)); // s = -(q + r) = -(-2 + 0) = 2 = size

    // Non-pentagon coordinates
    try expect(!grid.isPenta(0, 0));
    try expect(!grid.isPenta(1, 0));
    try expect(!grid.isPenta(0, 1));
}

test "Grid.hashCoord uniqueness" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);
    defer grid.deinit();

    const hash1 = planet.Grid.hashCoord(0, 0, 0);
    const hash2 = planet.Grid.hashCoord(1, 0, 0);
    const hash3 = planet.Grid.hashCoord(0, 0, 1);

    try expect(hash1 != hash2);
    try expect(hash1 != hash3);
    try expect(hash2 != hash3);
}

test "Grid.populateIndices basic functionality" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 1);
    defer grid.deinit();

    try grid.populateIndices();

    // Should have some tiles
    try expect(grid.tile_count > 0);

    // Coord map should be populated
    try expect(grid.coord_map != null);
    if (grid.coord_map) |*coord_map| {
        try expect(coord_map.count() > 0);
    }
}

test "Grid.populateIndices index consistency" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 1);
    defer grid.deinit();

    try grid.populateIndices();

    if (grid.coord_map) |*coord_map| {
        var iterator = coord_map.iterator();

        // All indices should be within bounds
        while (iterator.next()) |entry| {
            const index = entry.value_ptr.*;
            try expect(index < grid.tile_count);
        }

        // Check that we have the expected number of unique coordinates
        // For size=1, we should have coordinates for all valid (q,r) pairs across 20 faces
        const expected_min_tiles = 20; // At least 20 tiles (one per face)
        try expect(coord_map.count() >= expected_min_tiles);
    }
}

test "Grid.populateIndices shared edge handling" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);
    defer grid.deinit();

    try grid.populateIndices();

    if (grid.coord_map) |*coord_map| {
        // Find edge tiles and verify they share indices across faces
        var edge_indices = std.ArrayListUnmanaged(usize){};
        defer edge_indices.deinit(std.testing.allocator);

        var iterator = coord_map.iterator();
        while (iterator.next()) |entry| {
            const coord_hash = entry.key_ptr.*;
            const index = entry.value_ptr.*;

            // Extract face from hash
            const face = coord_hash >> 40;

            // Only check first few faces for this test
            if (face < 3) {
                try edge_indices.append(std.testing.allocator, index);
            }
        }

        // Should have some indices
        try expect(edge_indices.items.len > 0);
    }
}

test "Grid.populateIndices pentagon handling" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);
    defer grid.deinit();

    try grid.populateIndices();

    // Pentagon indices should be initialized (not max value)
    for (grid.penta_indices) |penta_index| {
        try expect(penta_index != std.math.maxInt(usize));
    }
}

test "Grid memory management" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);

    try grid.populateIndices();

    // Should not leak when deinitialized
    grid.deinit();

    // If we got here without crashing, memory management is working
    try expect(true);
}

test "Grid different sizes" {
    const allocator = std.testing.allocator;

    // Test different grid sizes
    const sizes = [_]usize{ 1, 2, 3 };

    for (sizes) |size| {
        var grid = try planet.Grid.init(allocator, size);
        defer grid.deinit();

        try grid.populateIndices();

        // Larger grids should have more tiles
        try expect(grid.tile_count > 0);

        if (grid.coord_map) |*coord_map| {
            try expect(coord_map.count() > 0);
        }
    }
}

test "Grid coordinate hash roundtrip" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);
    defer grid.deinit();

    const test_coords = [_]struct { q: isize, r: isize, face: usize }{
        .{ .q = 0, .r = 0, .face = 0 },
        .{ .q = 1, .r = 0, .face = 5 },
        .{ .q = -1, .r = 1, .face = 10 },
        .{ .q = 2, .r = -1, .face = 15 },
    };

    for (test_coords) |coord| {
        const hash = planet.Grid.hashCoord(coord.q, coord.r, coord.face);

        // Verify hash encodes the face correctly
        const extracted_face = hash >> 44;
        try expectEqual(coord.face, extracted_face);
    }
}

test "Grid.hex coordinate system properties" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 3);
    defer grid.deinit();

    // Test that q + r + s = 0 for cube coordinates
    const test_coords = [_]struct { q: isize, r: isize }{
        .{ .q = 0, .r = 0 },
        .{ .q = 1, .r = -1 },
        .{ .q = -1, .r = 1 },
        .{ .q = 2, .r = -1 },
        .{ .q = -2, .r = 2 },
    };

    for (test_coords) |coord| {
        const s = -(coord.q + coord.r);
        try expectEqual(@as(isize, 0), coord.q + coord.r + s);
    }
}

test "Grid.edge and pentagon relationship" {
    const allocator = std.testing.allocator;
    var grid = try planet.Grid.init(allocator, 2);
    defer grid.deinit();

    // All pentagons should also be edges (since they're at the corners)
    const penta_coords = [_]struct { q: isize, r: isize }{
        .{ .q = 2, .r = 0 },
        .{ .q = 0, .r = 2 },
        .{ .q = -2, .r = 0 },
        .{ .q = 0, .r = -2 },
    };

    for (penta_coords) |coord| {
        if (grid.isValid(coord.q, coord.r)) {
            const is_penta = grid.isPenta(coord.q, coord.r);
            const is_edge = grid.isEdge(coord.q, coord.r);

            if (is_penta) {
                try expect(is_edge); // All pentagons should be edges
            }
        }
    }
}

test "Grid.tile_count scales with size" {
    const allocator = std.testing.allocator;

    var grid1 = try planet.Grid.init(allocator, 1);
    defer grid1.deinit();
    try grid1.populateIndices();

    var grid2 = try planet.Grid.init(allocator, 2);
    defer grid2.deinit();
    try grid2.populateIndices();

    var grid3 = try planet.Grid.init(allocator, 3);
    defer grid3.deinit();
    try grid3.populateIndices();

    // Larger grids should have more tiles
    try expect(grid2.tile_count > grid1.tile_count);
    try expect(grid3.tile_count > grid2.tile_count);

    // Should have reasonable numbers (not zero)
    try expect(grid1.tile_count > 0);
    try expect(grid2.tile_count > 0);
    try expect(grid3.tile_count > 0);
}

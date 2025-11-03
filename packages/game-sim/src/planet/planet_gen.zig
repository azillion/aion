const std = @import("std");
const Grid = @import("grid.zig").Grid;
const SimplexNoise = @import("math").noise.SimplexNoise;

// Populate the grid's per-vertex elevations with FBM noise in [0,1].
pub fn generateBedrock(grid: *Grid, seed: u64) void {
    if (grid.vertices == null or grid.elevations == null) return;
    const verts = grid.vertices.?;
    const elevs = grid.elevations.?;

    const noise_gen = SimplexNoise.init(seed);

    var i: usize = 0;
    while (i < grid.tile_count) : (i += 1) {
        const base = i * 3;
        const x = @as(f64, @floatCast(verts[base + 0]));
        const y = @as(f64, @floatCast(verts[base + 1]));
        const z = @as(f64, @floatCast(verts[base + 2]));

        // Fractional Brownian Motion (multi-octave simplex)
        var total: f64 = 0.0;
        var frequency: f64 = 1.0;
        var amplitude: f64 = 0.5;
        const persistence: f64 = 0.5;
        const lacunarity: f64 = 2.0;

        var o: u32 = 0;
        while (o < 6) : (o += 1) {
            total += noise_gen.noise(x * frequency, y * frequency, z * frequency) * amplitude;
            frequency *= lacunarity;
            amplitude *= persistence;
        }

        // Normalize from roughly [-1,1] -> [0,1]
        const h = @as(f32, @floatCast(total * 0.5 + 0.5));
        elevs[i] = h;
    }
}

const expect = std.testing.expect;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;

test "planet_gen: bedrock elevations are [0,1] with variance" {
    var g = try Grid.init(std.testing.allocator, 3);
    defer g.deinit();

    var map = try g.populateIndices();
    defer map.deinit();
    try g.populateNeighbors(&map);
    try g.generateMesh();

    generateBedrock(&g, 12345);

    try expect(g.elevations != null);
    const elevs = g.elevations.?;
    var minv: f32 = 1e9;
    var maxv: f32 = -1e9;
    for (elevs) |h| {
        try expect(h >= 0.0 and h <= 1.0);
        if (h < minv) minv = h;
        if (h > maxv) maxv = h;
    }
    // Ensure we produced some non-constant structure
    try expect(maxv - minv > 0.01);
}

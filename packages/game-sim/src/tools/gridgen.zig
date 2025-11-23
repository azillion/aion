const std = @import("std");
const planet = @import("planet");
const Grid = planet.Grid;
const prebaked_builder = planet.prebaked_builder;
const prebaked_format = planet.prebaked_format;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    // Arbitrary resolution for testing the pipeline
    const R: u32 = 4;
    std.debug.print("GRIDGEN: Initializing icosahedral grid R={d}...\n", .{R});

    var grid = try Grid.init(alloc, R);
    defer grid.deinit();

    const planet_data = try prebaked_builder.build_prebaked(&grid, alloc);
    defer {
        alloc.free(planet_data.tiles);
        alloc.free(planet_data.vertices);
        alloc.free(planet_data.elevations);
        alloc.free(planet_data.indices);
        alloc.free(planet_data.edges);
    }

    std.debug.print(
        "GRIDGEN: Generated {d} tiles, {d} unique verts, {d} triangles.\n",
        .{ planet_data.tiles.len, planet_data.vertices.len / 3, planet_data.indices.len / 3 },
    );

    // Persist the prebaked data for runtime consumption.
    var name_buf: [32]u8 = undefined;
    const file_name = try std.fmt.bufPrint(&name_buf, "planet_R{d}.gridbin", .{R});
    var file = try std.fs.cwd().createFile(file_name, .{ .truncate = true });
    defer file.close();

    var header = prebaked_format.PlanetHeader{
        .R = R,
        .tile_count = @intCast(planet_data.tiles.len),
        .vertex_count = @intCast(planet_data.vertices.len / 3),
        .index_count = @intCast(planet_data.indices.len),
        .edge_index_count = @intCast(planet_data.edges.len),
    };
    try file.writeAll(std.mem.asBytes(&header));

    var i: usize = 0;
    while (i < planet_data.tiles.len) : (i += 1) {
        const disk = prebaked_format.tileToDisk(planet_data.tiles[i]);
        try file.writeAll(std.mem.asBytes(&disk));
    }
    try file.writeAll(std.mem.sliceAsBytes(planet_data.vertices));
    try file.writeAll(std.mem.sliceAsBytes(planet_data.elevations));
    try file.writeAll(std.mem.sliceAsBytes(planet_data.indices));
    try file.writeAll(std.mem.sliceAsBytes(planet_data.edges));

    std.debug.print("GRIDGEN: wrote {s} ({d} tiles)\n", .{ file_name, planet_data.tiles.len });

    // Sample a few tiles to prove coordinate mapping works
    if (planet_data.tiles.len > 0) {
        const t = planet_data.tiles[0];
        std.debug.print(
            "  Tile[0]: face={d} q={d} r={d} pos=({d:.3}, {d:.3}, {d:.3}) penta={any}\n",
            .{
                t.face,
                t.axial.q,
                t.axial.r,
                t.pos_sphere[0],
                t.pos_sphere[1],
                t.pos_sphere[2],
                t.is_pentagon,
            },
        );
        std.debug.print("    neighbors={any}\n", .{t.neighbors});
    }

    if (planet_data.tiles.len > 12) {
        const t = planet_data.tiles[12];
        std.debug.print(
            "  Tile[12]: face={d} q={d} r={d} pos=({d:.3}, {d:.3}, {d:.3}) penta={any}\n",
            .{
                t.face,
                t.axial.q,
                t.axial.r,
                t.pos_sphere[0],
                t.pos_sphere[1],
                t.pos_sphere[2],
                t.is_pentagon,
            },
        );
        std.debug.print("    neighbors={any}\n", .{t.neighbors});
    }
}

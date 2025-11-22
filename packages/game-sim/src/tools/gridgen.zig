const std = @import("std");
const planet = @import("planet");
const Grid = planet.Grid;
const prebaked_builder = planet.prebaked_builder;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const alloc = gpa.allocator();

    // Arbitrary resolution for testing the pipeline
    const R: u32 = 4;
    std.debug.print("GRIDGEN: Initializing icosahedral grid R={d}...\n", .{R});

    var grid = try Grid.init(alloc, R);
    defer grid.deinit();

    const tiles = try prebaked_builder.build_prebaked(&grid, alloc);
    defer alloc.free(tiles);

    std.debug.print("GRIDGEN: Generated {d} tiles.\n", .{tiles.len});

    // Sample a few tiles to prove coordinate mapping works
    if (tiles.len > 0) {
        const t = tiles[0];
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

    if (tiles.len > 12) {
        const t = tiles[12];
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

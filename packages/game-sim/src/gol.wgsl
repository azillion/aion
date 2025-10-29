@group(0) @binding(0)
var input_grid: texture_storage_2d<r32uint, read>;

@group(0) @binding(1)
var output_grid: texture_storage_2d<r32uint, write>;

fn wrap_coord(value: u32, delta: i32, max_value: u32) -> u32 {
    let vi = i32(value) + delta;
    // ensure positive modulo in range [0, max_value)
    return u32((vi + i32(max_value)) % i32(max_value));
}

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let dims: vec2<u32> = textureDimensions(input_grid);
    let width: u32 = dims.x;
    let height: u32 = dims.y;

    let x: u32 = gid.x;
    let y: u32 = gid.y;

    // Guard against over-dispatch just in case
    if (x >= width || y >= height) {
        return;
    }

    var live_neighbors: u32 = 0u;

    // Offsets for the 8 neighbors
    let offsets: array<vec2<i32>, 8> = array<vec2<i32>, 8>(
        vec2<i32>(-1, -1), vec2<i32>(0, -1), vec2<i32>(1, -1),
        vec2<i32>(-1,  0),                     vec2<i32>(1,  0),
        vec2<i32>(-1,  1), vec2<i32>(0,  1), vec2<i32>(1,  1)
    );

    for (var i: u32 = 0u; i < 8u; i = i + 1u) {
        let dx: i32 = offsets[i].x;
        let dy: i32 = offsets[i].y;
        let nx: u32 = wrap_coord(x, dx, width);
        let ny: u32 = wrap_coord(y, dy, height);
        let neighbor: u32 = textureLoad(input_grid, vec2<u32>(nx, ny)).r;
        live_neighbors = live_neighbors + (neighbor & 1u);
    }

    let current: u32 = textureLoad(input_grid, vec2<u32>(x, y)).r & 1u;

    // Conway's Game of Life rules
    var next: u32 = 0u;
    if (current == 1u) {
        if (live_neighbors == 2u || live_neighbors == 3u) {
            next = 1u;
        } else {
            next = 0u;
        }
    } else {
        if (live_neighbors == 3u) {
            next = 1u;
        } else {
            next = 0u;
        }
    }

    textureStore(output_grid, vec2<u32>(x, y), vec4<u32>(next, 0u, 0u, 0u));
}

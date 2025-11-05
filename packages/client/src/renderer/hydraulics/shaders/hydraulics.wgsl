#include "cubeSphere.wgsl"

const dt: f32 = 0.1;          // Simulation timestep
const flow_speed: f32 = 0.1;  // Flow coefficient

@group(0) @binding(0) var water_state_read: texture_2d_array<f32>;
@group(0) @binding(1) var water_state_write: texture_storage_2d_array<rgba16float, write>;
@group(0) @binding(2) var terrain_height: texture_2d_array<f32>;

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let dims = textureDimensions(water_state_read);
    if (id.x >= dims.x || id.y >= dims.y) {
        return;
    }

    let coords = vec2i(id.xy);
    let face_index = id.z;

    // --- Read previous state ---
    let water_h = textureLoad(water_state_read, coords, i32(face_index), 0).r;
    let terrain_h = textureLoad(terrain_height, coords, i32(face_index), 0).r;
    let H_me = terrain_h + water_h;

    var net_flow = 0.0;

    // --- Seam-aware neighbor flow ---
    let neighbors = array<vec2i, 4>(vec2i(1,0), vec2i(-1,0), vec2i(0,1), vec2i(0,-1));
    for (var i = 0; i < 4; i = i + 1) {
        let nbr = getNeighbor(coords, face_index, neighbors[i], dims.x);
        let n_coords = vec2i(nbr.xy);
        let n_face = i32(nbr.z);
        let water_h_n = textureLoad(water_state_read, n_coords, n_face, 0).r;
        let terrain_h_n = textureLoad(terrain_height, n_coords, n_face, 0).r;
        let H_n = terrain_h_n + water_h_n;
        net_flow += max(0.0, H_n - H_me) * flow_speed;
        net_flow -= max(0.0, H_me - H_n) * flow_speed;
    }
    var new_water_h = max(0.0, water_h + net_flow * dt);

    // Light rain input
    let rain_rate: f32 = 0.00001;
    new_water_h += rain_rate * dt;

    // --- Write new state ---
    textureStore(water_state_write, coords, i32(face_index), vec4f(new_water_h, 0.0, 0.0, 1.0));
}



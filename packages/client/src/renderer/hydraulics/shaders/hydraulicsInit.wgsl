@group(0) @binding(0) var water_state_write: texture_storage_2d_array<rgba16float, write>;
@group(0) @binding(1) var terrain_height: texture_2d_array<f32>;

struct PlanetParams {
    radius: f32,
    sea_level: f32,
};
@group(0) @binding(2) var<uniform> planet: PlanetParams;

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let dims = textureDimensions(water_state_write);
    if (id.x >= dims.x || id.y >= dims.y) { return; }
    let coords = vec2i(id.xy);
    let face = i32(id.z);
    let terrain_h : f32 = textureLoad(terrain_height, coords, face, 0).r;
    let initial_water = max(0.0, planet.sea_level - terrain_h);
    textureStore(water_state_write, coords, face, vec4f(initial_water, 0.0, 0.0, 1.0));
}



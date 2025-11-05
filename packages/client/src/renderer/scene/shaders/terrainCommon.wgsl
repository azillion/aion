// Common terrain-related uniforms and types shared between terrain and shading.
// This must match the layout of the terrain_params vec4 in the Sphere struct.
struct TerrainUniforms {
    base_radius: f32,
    sea_level: f32,
    max_height: f32,
    seed: f32,
};



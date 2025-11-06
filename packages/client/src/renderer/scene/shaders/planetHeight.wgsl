// Minimal terrain height utilities for baking.
// Depends on snoise (noise.wgsl) and sampleCoarseGrid (coarseGrid.wgsl) being included by the entry shader.

fn warp_bake(p: vec3<f32>, seed: f32) -> vec3<f32> {
    let q = vec3<f32>(
        snoise(p + vec3<f32>(0.0, 0.0, seed)),
        snoise(p + vec3<f32>(5.2, 1.3, seed)),
        snoise(p + vec3<f32>(2.1, 3.4, seed))
    );
    return p + 0.8 * q;
}

fn fbm_bake(p: vec3<f32>, seed: f32, octaves: i32) -> f32 {
    var total: f32 = 0.0;
    var frequency: f32 = 4.0;
    var amplitude: f32 = 0.5;
    let persistence: f32 = 0.5;
    let lacunarity: f32 = 2.0;

    for (var i = 0; i < octaves; i = i + 1) {
        total = total + snoise(p * frequency + seed) * amplitude;
        frequency = frequency * lacunarity;
        amplitude = amplitude * persistence;
    }
    return total;
}

// Calculates signed terrain height (km) based on coarse grid + FBM detail with LOD.
fn h_noise(dir: vec3f, params: TerrainUniforms, dist_to_surface: f32, base_radius: f32, camera: CameraUniforms, scene: SceneUniforms) -> f32 {
    // 1) Base elevation from coarse grid [0,1] - returns a single f32
    let base_elevation_norm = sampleCoarseGrid(dir);

    // 2) Detail noise with LOD
    // Smoothly transitions from 1 octave (far) to 8 octaves (very close)
    let lod_f = clamp(1.0 - (log(dist_to_surface / base_radius) / log(10000.0)), 0.0, 1.0);
    let num_octaves = i32(1.0 + lod_f * 7.0);
    let detail_noise = fbm_bake(dir, params.seed, num_octaves);

    // 3) Combine
    let max_height_km = params.max_height * params.base_radius;
    let base_height_km = (base_elevation_norm - 0.5) * max_height_km;
    let detail_height_km = detail_noise * (max_height_km * 0.1) * lod_f; // fade detail at distance
    return base_height_km + detail_height_km;
}

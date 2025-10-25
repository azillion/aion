// This must match the layout of the terrain_params vec4 in the Sphere struct.
struct TerrainUniforms {
    base_radius: f32,
    sea_level: f32,
    max_height: f32,
    seed: f32,
};

// Domain warping function for more interesting shapes.
fn warp(p: vec3<f32>, seed: f32) -> vec3<f32> {
    let q = vec3<f32>(
        snoise(p + vec3<f32>(0.0, 0.0, seed)),
        snoise(p + vec3<f32>(5.2, 1.3, seed)),
        snoise(p + vec3<f32>(2.1, 3.4, seed))
    );
    return p + 0.8 * q;
}

// Calculates signed terrain height (km) using Fractal Brownian Motion (FBM).
fn h_noise(dir: vec3<f32>, params: TerrainUniforms, dist_to_surface: f32, base_radius: f32, camera: CameraUniforms) -> f32 {
    var h = 0.0;
    var a = 1.0;
    var f = 2.5;
    let octaves = 6; // Fixed octaves for now.

    var total_amplitude = 0.0;
    for(var i = 0; i < octaves; i = i + 1) {
        // Derived LOD: stop when projected feature size < 1px (with quality factor)
        let quality_factor = 2.0;
        if (dist_to_surface > (base_radius / f) * (camera.projection_constants.x / quality_factor)) { break; }
        let p = dir * f;
        h = h + a * snoise(p);
        total_amplitude = total_amplitude + a;
        a = a * 0.5;
        f = f * 2.0;
    }

    if (total_amplitude == 0.0) { return 0.0; }
    let normalized_h = h / total_amplitude;
    // Map [-1,1] -> [0,1], then scale by base_radius and proportional max_height scalar
    let height_in_km = (normalized_h * 0.5 + 0.5) * params.base_radius * params.max_height;
    return height_in_km;
}

// Calculates a more accurate normal for a displaced heightfield via central differencing.
// p is a point on the analytic sphere surface in world space.
// planet_center ensures we compute directions relative to the planet, not world origin.
fn get_terrain_normal_from_heightfield(
    p_on_sphere: vec3<f32>,
    planet_center: vec3<f32>,
    analytic_normal: vec3<f32>,
    params: TerrainUniforms,
    dist_to_surface: f32,
    base_radius: f32,
    camera: CameraUniforms
) -> vec3<f32> {
    // Adaptive epsilon proportional to camera distance to avoid catastrophic cancellation
    let eps = max(0.001, dist_to_surface * 0.0001);

    var up_vec = vec3<f32>(0.0, 1.0, 0.0);
    if (abs(dot(analytic_normal, up_vec)) > 0.999) { up_vec = vec3<f32>(1.0, 0.0, 0.0); }
    let tangent = normalize(cross(analytic_normal, up_vec));
    let bitangent = normalize(cross(analytic_normal, tangent));

    let h0 = h_noise(normalize(p_on_sphere - planet_center), params, dist_to_surface, base_radius, camera);
    let h1 = h_noise(normalize(p_on_sphere + tangent * eps - planet_center), params, dist_to_surface, base_radius, camera);
    let h2 = h_noise(normalize(p_on_sphere + bitangent * eps - planet_center), params, dist_to_surface, base_radius, camera);

    let final_normal = normalize(
        analytic_normal - tangent * (h1 - h0) / eps - bitangent * (h2 - h0) / eps
    );
    return final_normal;
}


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
    var f = 4.0;
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
    camera: CameraUniforms,
    scene: SceneUniforms
) -> vec3<f32> {
    // Adaptive epsilon proportional to camera distance to avoid catastrophic cancellation
    let eps = max(0.001, dist_to_surface * 0.00001);
    let tier_scale = scene.tier_scale_and_pad.x;

    var up_vec = vec3<f32>(0.0, 1.0, 0.0);
    if (abs(dot(analytic_normal, up_vec)) > 0.999) { up_vec = vec3<f32>(1.0, 0.0, 0.0); }
    let tangent = normalize(cross(analytic_normal, up_vec));
    let bitangent = normalize(cross(analytic_normal, tangent));

    let h0 = h_noise(normalize(p_on_sphere - planet_center), params, dist_to_surface, base_radius, camera);
    let h1 = h_noise(normalize(p_on_sphere + tangent * eps - planet_center), params, dist_to_surface + eps, base_radius, camera);
    let h2 = h_noise(normalize(p_on_sphere + bitangent * eps - planet_center), params, dist_to_surface + eps, base_radius, camera);

    let final_normal = normalize(
        analytic_normal - (tangent * (h1 - h0) / eps + bitangent * (h2 - h0) / eps) * (1.0 / tier_scale)
    );
    return final_normal;
}

fn get_sdf_normal(
    p: vec3<f32>,
    planet_center: vec3<f32>,
    params: TerrainUniforms,
    camera: CameraUniforms,
    scene: SceneUniforms
) -> vec3<f32> {
    let eps = 0.01;
    let dx = vec3<f32>(eps, 0.0, 0.0);
    let dy = vec3<f32>(0.0, eps, 0.0);
    let dz = vec3<f32>(0.0, 0.0, eps);

    let p_local = p - planet_center;
    let dist_to_p = length(p);

    let nx = dWorld(p_local + dx, params, dist_to_p, camera, scene) - dWorld(p_local - dx, params, dist_to_p, camera, scene);
    let ny = dWorld(p_local + dy, params, dist_to_p, camera, scene) - dWorld(p_local - dy, params, dist_to_p, camera, scene);
    let nz = dWorld(p_local + dz, params, dist_to_p, camera, scene) - dWorld(p_local - dz, params, dist_to_p, camera, scene);

    return normalize(vec3<f32>(nx, ny, nz));
}

fn dWorld(
    p_local: vec3<f32>,
    params: TerrainUniforms,
    dist_marched: f32,
    camera: CameraUniforms,
    scene: SceneUniforms
) -> f32 {
    let dir = normalize(p_local);
    let tier_scale = scene.tier_scale_and_pad.x;
    let real_dist_for_lod = dist_marched * tier_scale;
    let h = h_noise(dir, params, real_dist_for_lod, params.base_radius, camera);
    let h_scaled = h / tier_scale;
    let R_scaled = params.base_radius / tier_scale;
    let d_terrain = length(p_local) - (R_scaled + h_scaled);
    let water_radius_scaled = (params.base_radius + params.sea_level) / tier_scale;
    let d_ocean = length(p_local) - water_radius_scaled;
    return min(d_terrain, d_ocean);
}

fn ray_march(
    ray: Ray,
    max_dist: f32,
    planet_center: vec3<f32>,
    params: TerrainUniforms,
    camera: CameraUniforms,
    scene: SceneUniforms
) -> HitRecord {
    var rec: HitRecord; rec.hit = false;
    var t = 0.0;

    for (var i = 0; i < 128; i = i + 1) {
        let p = rayAt(ray, t);
        let p_local = p - planet_center;
        let d = dWorld(p_local, params, t, camera, scene);

        // Adaptive hit threshold proportional to distance marched
        let threshold = max(0.001, 0.001 * t);
        if (d < threshold) {
            rec.hit = true;
            rec.t = t;
            rec.p = p;
            let outward_normal = get_sdf_normal(p, planet_center, params, camera, scene);
            rec.front_face = dot(ray.direction, outward_normal) < 0.0;
            rec.normal = outward_normal;
            return rec;
        }

        t = t + max(0.001, d * 0.8);
        if (t > max_dist) { break; }
    }

    return rec;
}


// --- Phase 3: Shading & Materials ---

struct Material {
    albedo: vec3<f32>,
};

// Simple procedural texture returning grayscale variation
fn tex3(p: vec3<f32>) -> vec3<f32> {
    let scale = 0.05;
    let n = snoise(p * scale) * 0.5 + 0.5;
    return vec3<f32>(n);
}

// 2D variant used by triplanar projections
fn tex2(p: vec2<f32>) -> vec3<f32> {
    let scale = 0.05;
    let n = snoise(vec3<f32>(p * scale, 0.0)) * 0.5 + 0.5;
    return vec3<f32>(n);
}

// Triplanar mapping with normalized blend weights
fn tex_triplanar(p: vec3<f32>, normal: vec3<f32>) -> vec3<f32> {
    let weights = abs(normal);
    let denom = max(1e-5, weights.x + weights.y + weights.z);
    let blend = weights / denom;

    let x_proj = tex2(p.yz);
    let y_proj = tex2(p.xz);
    let z_proj = tex2(p.xy);

    return x_proj * blend.x + y_proj * blend.y + z_proj * blend.z;
}

fn get_material(p_local: vec3<f32>, normal: vec3<f32>, params: TerrainUniforms, h_scaled: f32) -> Material {
    let altitude = h_scaled;
    let slope = 1.0 - dot(normal, normalize(p_local));

    // Base material colors modulated by triplanar texture
    let texval = tex_triplanar(p_local, normal);
    let rock_color = vec3<f32>(0.5, 0.45, 0.4) * texval;
    let ground_color = vec3<f32>(0.3, 0.4, 0.15) * texval;
    let snow_color = vec3<f32>(0.9, 0.9, 0.95) * texval;

    // Blends
    let rock_amount = smoothstep(0.4, 0.6, slope);
    let snow_start = params.max_height * params.base_radius * 0.6;
    let snow_end = params.max_height * params.base_radius * 0.8;
    let snow_amount = smoothstep(snow_start, snow_end, altitude);

    var albedo = mix(ground_color, rock_color, rock_amount);
    albedo = mix(albedo, snow_color, snow_amount);
    return Material(albedo);
}

fn get_ocean_wave_normal(p_local: vec3<f32>, base_normal: vec3<f32>) -> vec3<f32> {
    let wave_freq = 50.0;
    let wave_amp = 0.01;
    let wave_noise = snoise(p_local * wave_freq);
    return normalize(base_normal + wave_amp * vec3<f32>(wave_noise));
}

fn get_ocean_material(p_local: vec3<f32>) -> Material {
    // Base, unlit ocean color
    let water_albedo = vec3<f32>(0.05, 0.15, 0.2);
    return Material(water_albedo);
}


const ATMOSPHERE_RADIUS_SCALE: f32 = 1.08;
const PI: f32 = 3.1415926535;

struct AtmosphereParams {
    planet_radius: f32,
    atmosphere_radius: f32,
    h_rayleigh: f32,
    h_mie: f32,
    beta_rayleigh: vec3<f32>,
    beta_mie: vec3<f32>,
    g_mie: f32,
};

fn get_earth_atmosphere() -> AtmosphereParams {
    return AtmosphereParams(
        6371.0,
        6371.0 * ATMOSPHERE_RADIUS_SCALE,
        8.5,
        1.2,
        vec3<f32>(5.802e-3, 9.558e-3, 33.1e-3),
        vec3<f32>(3.996e-3),
        0.8
    );
}

fn ray_sphere_intersect(ray_origin: vec3<f32>, ray_dir: vec3<f32>, sphere_radius: f32) -> vec2<f32> {
    let b = dot(ray_origin, ray_dir);
    let c = dot(ray_origin, ray_origin) - sphere_radius * sphere_radius;
    var d = b * b - c;
    if (d < 0.0) { return vec2<f32>(-1.0); }
    d = sqrt(d);
    return vec2<f32>(-b - d, -b + d);
}

struct AtmosphereOutput {
    in_scattering: vec3<f32>,
    transmittance: vec3<f32>,
    alpha: f32,
};

fn get_optical_depth(
    p_sample: vec3<f32>,
    light_dir: vec3<f32>,
    params: AtmosphereParams,
    tier_scale: f32,
    planet_radius_scaled: f32,
    atmosphere_radius_scaled: f32
) -> vec3<f32> {
    let t_light_atmos = ray_sphere_intersect(p_sample, light_dir, atmosphere_radius_scaled);
    let path_length_to_atmos_edge = t_light_atmos.y;
    if (path_length_to_atmos_edge <= 0.0) { return vec3<f32>(1e9); }
    let num_light_samples = 8;
    let step_size_scaled = path_length_to_atmos_edge / f32(num_light_samples);
    var optical_depth = vec3<f32>(0.0);
    for (var j = 0; j < num_light_samples; j = j + 1) {
        let t_j = (f32(j) + 0.5) * step_size_scaled;
        let p_light_sample = p_sample + light_dir * t_j;
        if (length(p_light_sample) < planet_radius_scaled) {
            optical_depth += vec3<f32>(1e9);
            continue;
        }
        let h_light = max(0.0, length(p_light_sample) * tier_scale - params.planet_radius);
        let density_rayleigh = exp(-h_light / params.h_rayleigh);
        let density_mie = exp(-h_light / params.h_mie);
        let step_size_world = step_size_scaled * tier_scale;
        optical_depth += (params.beta_rayleigh * density_rayleigh + params.beta_mie * density_mie) * step_size_world;
    }
    return optical_depth;
}

// fn get_sky_color(ray: Ray, planet_center_world: vec3<f32>, planet_radius_scaled: f32, light_dir: vec3<f32>, scene: SceneUniforms, max_distance: f32) -> AtmosphereOutput {
//     // The camera is at (0,0,0), so this vector is FROM the camera TO the planet center.
//     let ray_origin_local = -planet_center_world;
//     let dist_to_center = length(ray_origin_local);

//     // Use the provided scaled planet radius.
//     let atmosphere_radius_scaled = planet_radius_scaled * ATMOSPHERE_RADIUS_SCALE;

//     // --- THE TEST ---
//     // If the GPU thinks the camera is INSIDE the atmosphere, the screen will be GREEN.
//     // If the GPU thinks the camera is OUTSIDE the atmosphere, the screen will be RED.
//     if (dist_to_center < atmosphere_radius_scaled) {
//         return AtmosphereOutput(vec3<f32>(0.0, 1.0, 0.0), vec3<f32>(0.0), 1.0); // GREEN = INSIDE
//     } else {
//         return AtmosphereOutput(vec3<f32>(1.0, 0.0, 0.0), vec3<f32>(0.0), 1.0); // RED = OUTSIDE
//     }
// }

fn get_sky_color(ray: Ray, planet_center_world: vec3<f32>, planet_radius_scaled: f32, light_dir: vec3<f32>, scene: SceneUniforms, max_distance: f32) -> AtmosphereOutput {
    let params = get_earth_atmosphere();
    let tier_scale = scene.tier_scale_and_pad.x;
    let atmosphere_radius_scaled = planet_radius_scaled * ATMOSPHERE_RADIUS_SCALE;
    let ray_origin_local = -planet_center_world;
    let t = ray_sphere_intersect(ray_origin_local, ray.direction, atmosphere_radius_scaled);
    if (t.y < 0.0) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0), 0.0); }
    let start_offset = max(0.0, t.x);
    let end_t = select(t.y, min(t.y, max_distance), max_distance > 0.0);
    let dist_in_atmosphere = end_t - start_offset;
    if (dist_in_atmosphere <= 0.0) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0), 0.0); }
    let num_samples = 16;
    let step_size_scaled = dist_in_atmosphere / f32(num_samples);
    let step_size_world = step_size_scaled * tier_scale;
    var transmittance = vec3<f32>(1.0);
    var scattered_light = vec3<f32>(0.0);
    for (var i = 0; i < num_samples; i = i + 1) {
        let t_sample = start_offset + (f32(i) + 0.5) * step_size_scaled;
        let p_sample = ray_origin_local + ray.direction * t_sample;
        let h = max(0.0, length(p_sample) * tier_scale - params.planet_radius);
        let density_rayleigh = exp(-h / params.h_rayleigh);
        let density_mie = exp(-h / params.h_mie);
        let scattering_coefficients = params.beta_rayleigh * density_rayleigh + params.beta_mie * density_mie;
        let optical_depth_to_sun = get_optical_depth(p_sample, light_dir, params, tier_scale, planet_radius_scaled, atmosphere_radius_scaled);
        let sun_transmittance = exp(-optical_depth_to_sun);
        let mu = dot(ray.direction, light_dir);
        let phase_rayleigh = 3.0 / (16.0 * PI) * (1.0 + mu * mu);
        let g = params.g_mie; let g2 = g * g;
        let phase_mie = (3.0 * (1.0 - g*g)) / (8.0 * PI * (2.0 + g*g)) * (1.0 + mu*mu) / max(1e-6, pow(1.0 + g*g - 2.0*g*mu, 1.5));
        let in_scatter = (params.beta_rayleigh * density_rayleigh * phase_rayleigh) + (params.beta_mie * density_mie * phase_mie);
        scattered_light += transmittance * in_scatter * sun_transmittance * step_size_world;
        transmittance *= exp(-scattering_coefficients * step_size_world);
    }
    let final_color = scattered_light * 20.0 * scene.dominant_light_color_and_debug.xyz;
    let brightness = dot(final_color, vec3<f32>(0.2126, 0.7152, 0.0722));
    let alpha = 1.0 - exp(-brightness * 2.0);
    return AtmosphereOutput(final_color, transmittance, alpha);
}
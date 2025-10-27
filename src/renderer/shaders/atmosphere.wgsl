const ATMOSPHERE_RADIUS_SCALE: f32 = 1.025;
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
        // Coefficients expressed per-kilometer to match world-space units (km)
        // Tuned lower to increase transmittance and reduce overly dark appearance
        vec3<f32>(2.0e-3, 4.5e-3, 11.0e-3),
        vec3<f32>(0.2e-3),
        0.76
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
    p_start_local: vec3<f32>,     // Starting point in local tier-space (relative to planet center)
    light_dir: vec3<f32>,         // World-space direction TO the sun
    params: AtmosphereParams,
    tier_scale: f32,
    planet_radius_local: f32,     // Planet's solid radius in local tier-space
    atmosphere_radius_local: f32  // Atmosphere's outer radius in local tier-space
) -> vec3<f32> {
    // 1) Intersect light ray with atmosphere boundary to get maximum path length
    let t_to_atmos_edge = ray_sphere_intersect(p_start_local, light_dir, atmosphere_radius_local).y;
    if (t_to_atmos_edge <= 0.0) { return vec3<f32>(0.0); }

    // 2) March along path; check shadowing per-sample
    let num_light_samples = 8;
    let step_size_local = t_to_atmos_edge / f32(num_light_samples);
    var optical_depth = vec3<f32>(0.0);

    for (var j = 0; j < num_light_samples; j = j + 1) {
        let p_sample_local = p_start_local + light_dir * (f32(j) + 0.5) * step_size_local;

        // Shadowing: inside solid planet blocks sun entirely (use geometric radius)
        if (length(p_sample_local) < planet_radius_local) { return vec3<f32>(1e9); }

        // Accumulate density using reconstructed physical radius for altitude based on the physical planet radius
        let altitude_km = max(0.0, length(p_sample_local) * tier_scale - params.planet_radius);
        if (altitude_km > 200.0) { continue; }

        let density_rayleigh = exp(-altitude_km / params.h_rayleigh);
        let density_mie = exp(-altitude_km / params.h_mie);
        let step_size_world = step_size_local * tier_scale;
        optical_depth += (params.beta_rayleigh * density_rayleigh + params.beta_mie * density_mie) * step_size_world;
    }
    return optical_depth;
}

fn get_sky_color(ray: Ray, planet_center_world: vec3<f32>, planet_radius_local: f32, light_dir: vec3<f32>, scene: SceneUniforms, max_distance: f32) -> AtmosphereOutput {
    let params = get_earth_atmosphere();
    let tier_scale = scene.tier_scale_and_pad.x;
    // Derive atmosphere radius in local tier-space from the local geometric planet radius
    let atmosphere_radius_local = planet_radius_local * ATMOSPHERE_RADIUS_SCALE;
    let ray_origin_local = ray.origin - planet_center_world;
    let t = ray_sphere_intersect(ray_origin_local, ray.direction, atmosphere_radius_local);
    
    // If the ray doesn't intersect the atmosphere at all, return nothing. This fixes the "horns".
    if (t.y <= 0.0) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0), 0.0); }
    
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

        // Don't accumulate light from under the surface
        if (length(p_sample) < planet_radius_local) { continue; }

        // Use reconstructed physical altitude above the true planet radius
        let h = max(0.0, length(p_sample) * tier_scale - params.planet_radius);
        let density_rayleigh = exp(-h / params.h_rayleigh);
        let density_mie = exp(-h / params.h_mie);
        let scattering_coefficients = params.beta_rayleigh * density_rayleigh + params.beta_mie * density_mie;
        let optical_depth_to_sun = get_optical_depth(p_sample, light_dir, params, tier_scale, planet_radius_local, atmosphere_radius_local);
        let sun_transmittance = exp(-optical_depth_to_sun);
        // Smoothly suppress scattering on the night side (slightly softened terminator)
        let day_night_falloff = smoothstep(-0.2, 0.05, dot(normalize(p_sample), light_dir));
        let mu = dot(ray.direction, light_dir);
        let phase_rayleigh = 3.0 / (16.0 * PI) * (1.0 + mu * mu);
        let g = params.g_mie; let g2 = g * g;
        let phase_mie = (3.0 * (1.0 - g*g)) / (8.0 * PI * (2.0 + g*g)) * (1.0 + mu*mu) / max(1e-6, pow(1.0 + g*g - 2.0*g*mu, 1.5));
        let in_scatter = (params.beta_rayleigh * density_rayleigh * phase_rayleigh) + (params.beta_mie * density_mie * phase_mie);
        scattered_light += transmittance * in_scatter * sun_transmittance * day_night_falloff * step_size_world;
        transmittance *= exp(-scattering_coefficients * step_size_world);
    }

    // // --- TEMPORARY DEBUG VIEW ---
    // // Visualize the sun color through the atmosphere at the midpoint of the segment.
    // let t_mid = start_offset + dist_in_atmosphere * 0.5;
    // let debug_sample_point = ray_origin_local + ray.direction * t_mid;
    // let optical_depth_to_sun_debug = get_optical_depth(debug_sample_point, light_dir, params, tier_scale, planet_radius_scaled, atmosphere_radius_scaled);
    // let sun_transmittance_debug = exp(-optical_depth_to_sun_debug);
    // return AtmosphereOutput(sun_transmittance_debug * scene.dominant_light_color_and_debug.xyz, vec3<f32>(1.0), 1.0);
    // // --- END DEBUG VIEW ---

    let final_color = scattered_light * 45.0 * scene.dominant_light_color_and_debug.xyz;
    let brightness = dot(final_color, vec3<f32>(0.2126, 0.7152, 0.0722));
    let alpha = 1.0 - exp(-brightness * 2.0);
    return AtmosphereOutput(final_color, transmittance, alpha);
}
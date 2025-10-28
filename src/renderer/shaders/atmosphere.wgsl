const ATMOSPHERE_RADIUS_SCALE: f32 = 1.025;

struct AtmosphereParams {
    planet_radius: f32,
    atmosphere_radius: f32,
    h_rayleigh: f32, // Rayleigh scale height (km)
    h_mie: f32,      // Mie scale height (km)
    beta_rayleigh: vec3<f32>, // Scattering coefficients at sea level
    beta_mie: vec3<f32>,
    g_mie: f32,      // Mie scattering asymmetry
};

// Earth-like parameters in world units (kilometers)
fn get_earth_atmosphere() -> AtmosphereParams {
    return AtmosphereParams(
        6371.0, // planet_radius
        6371.0 * ATMOSPHERE_RADIUS_SCALE, // atmosphere_radius
        8.5,    // h_rayleigh
        1.2,    // h_mie
        vec3<f32>(5.8e-3, 13.5e-3, 33.1e-3), // beta_rayleigh (tuned for visual appeal)
        vec3<f32>(2.1e-3), // beta_mie
        0.76    // g_mie
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

// Calculates optical depth for a ray segment through the atmosphere.
fn get_optical_depth(
    p_start_local: vec3<f32>,
    light_dir: vec3<f32>,
    params: AtmosphereParams,
    tier_scale: f32,
    planet_radius_local: f32,
    atmosphere_radius_local: f32,
    physical_radius_world: f32
) -> vec3<f32> {
    let t_to_atmos_edge = ray_sphere_intersect(p_start_local, light_dir, atmosphere_radius_local).y;
    if (t_to_atmos_edge <= 0.0) { return vec3<f32>(0.0); }

    let num_light_samples = 8;
    let step_size_local = t_to_atmos_edge / f32(num_light_samples);
    var optical_depth = vec3<f32>(0.0);
    let R = physical_radius_world;
    let R2 = R * R;

    for (var j = 0; j < num_light_samples; j = j + 1) {
        let p_sample_local = p_start_local + light_dir * (f32(j) + 0.5) * step_size_local;
        if (length(p_sample_local) < planet_radius_local) { return vec3<f32>(1e9); } // Shadowed by planet

        // Numerically stable altitude calculation for points far from the origin
        let d2 = dot(p_sample_local, p_sample_local) * tier_scale * tier_scale;
        let h = (d2 - R2) / (2.0 * R);
        
        if (h < 0.0 || h > 200.0) { continue; } // Sample is outside effective atmosphere
        
        let density_r = exp(-h / params.h_rayleigh);
        let density_m = exp(-h / params.h_mie);
        optical_depth += (params.beta_rayleigh * density_r + params.beta_mie * density_m) * (step_size_local * tier_scale);
    }
    return optical_depth;
}

struct AtmosphereOutput {
    in_scattering: vec3<f32>,
    transmittance: vec3<f32>,
    alpha: f32,
};

// Calculates scattered light and transmittance for a camera ray.
fn get_sky_color(
    ray: Ray, // in local space, relative to planet center
    planet_radius_local: f32,
    physical_radius_world: f32,
    light_dir_to_sun: vec3<f32>,
    light_color: vec3<f32>,
    scene: SceneUniforms,
    max_dist_local: f32
) -> AtmosphereOutput {
    let params = get_earth_atmosphere();
    let tier_scale = scene.tier_scale_and_pad.x;
    let atmosphere_radius_local = planet_radius_local * ATMOSPHERE_RADIUS_SCALE;
    let R = physical_radius_world;
    let R2 = R * R;
    
    let t_atmos = ray_sphere_intersect(ray.origin, ray.direction, atmosphere_radius_local);
    if (t_atmos.y < 0.0) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0), 0.0); }

    let start_t = max(0.0, t_atmos.x);
    let end_t = select(t_atmos.y, min(t_atmos.y, max_dist_local), max_dist_local < INFINITY);
    let dist_in_atmos = end_t - start_t;
    if (dist_in_atmos <= 0.0) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0), 0.0); }
    
    let num_samples = 16;
    let step_size_local = dist_in_atmos / f32(num_samples);
    var transmittance = vec3<f32>(1.0);
    var scattered_light = vec3<f32>(0.0);

    for (var i = 0; i < num_samples; i = i + 1) {
        let p_sample_local = ray.origin + (start_t + (f32(i) + 0.5) * step_size_local) * ray.direction;
        if (length(p_sample_local) < planet_radius_local) { continue; }

        let d2 = dot(p_sample_local, p_sample_local) * tier_scale * tier_scale;
        let h = (d2 - R2) / (2.0 * R);
        if (h < 0.0 || h > 200.0) { continue; }

        let density_r = exp(-h / params.h_rayleigh);
        let density_m = exp(-h / params.h_mie);
        let scattering_coeffs = params.beta_rayleigh * density_r + params.beta_mie * density_m;
        
        let optical_depth_to_sun = get_optical_depth(p_sample_local, light_dir_to_sun, params, tier_scale, planet_radius_local, atmosphere_radius_local, physical_radius_world);
        let sun_transmittance = exp(-optical_depth_to_sun);
        
        let mu = dot(ray.direction, light_dir_to_sun);
        let phase_rayleigh = 3.0 / (16.0 * PI) * (1.0 + mu * mu);
        let g = params.g_mie;
        let g2 = g * g;
        let phase_mie = (3.0 * (1.0 - g2)) / (8.0 * PI * (2.0 + g2)) * (1.0 + mu * mu) / pow(1.0 + g2 - 2.0 * g * mu, 1.5);
        
        let in_scatter = (params.beta_rayleigh * density_r * phase_rayleigh) + (params.beta_mie * density_m * phase_mie);
        
        scattered_light += transmittance * in_scatter * sun_transmittance * (step_size_local * tier_scale);
        transmittance *= exp(-scattering_coeffs * (step_size_local * tier_scale));
    }

    let final_color = scattered_light * light_color;
    let brightness = dot(final_color, vec3<f32>(0.2126, 0.7152, 0.0722));
    let alpha = 1.0 - exp(-brightness * 2.0); // Simple alpha for blending skies
    return AtmosphereOutput(final_color, transmittance, alpha);
}

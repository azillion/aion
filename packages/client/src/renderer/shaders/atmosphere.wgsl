const ATMOSPHERE_RADIUS_SCALE: f32 = 1.025;
const NUM_IN_SCATTER_SAMPLES: u32 = 16;

struct AtmosphereParams {
    planet_radius: f32,
    atmosphere_radius: f32,
    
    // Rayleigh
    beta_rayleigh: vec3<f32>,
    h_rayleigh: f32,
    
    // Mie
    beta_mie: vec3<f32>,
    h_mie: f32,
    g_mie: f32,
};

fn get_earth_atmosphere() -> AtmosphereParams {
    return AtmosphereParams(
        6371.0, 6371.0 * ATMOSPHERE_RADIUS_SCALE,
        vec3<f32>(5.8e-3, 13.5e-3, 33.1e-3), 8.5,
        vec3<f32>(0.4e-3), 1.2,
        0.76
    );
}

fn ray_sphere_intersect(r0: vec3<f32>, rd: vec3<f32>, sr: f32) -> vec2<f32> {
    let b = dot(rd, r0);
    let c = dot(r0, r0) - (sr * sr);
    let d = b * b - c;
    if (d < 0.0) { return vec2<f32>(1e5, -1e5); }
    return vec2<f32>(-b - sqrt(d), -b + sqrt(d));
}

struct AtmosphereOutput { in_scattering: vec3<f32>, transmittance: vec3<f32>, alpha: f32 };

fn calculate_optical_length_to_space(
    p_start: vec3<f32>,
    ray_dir: vec3<f32>,
    planet_radius: f32,
    atmosphere_radius: f32,
    h_rayleigh: f32,
    h_mie: f32,
    tier_scale: f32
) -> vec2<f32> {
    let t_atmos_exit = ray_sphere_intersect(p_start, ray_dir, atmosphere_radius).y;
    if (t_atmos_exit <= 0.0) {
        return vec2<f32>(0.0);
    }

    let t_planet_intersect = ray_sphere_intersect(p_start, ray_dir, planet_radius).x;
    if (t_planet_intersect > 0.0 && t_planet_intersect < t_atmos_exit) {
        return vec2<f32>(1e9); // Return a huge optical length for full occlusion.
    }
    
    let path_length = t_atmos_exit;
    let p_mid = p_start + ray_dir * (path_length * 0.5);

    let h_mid_local = length(p_mid) - planet_radius;
    var h_mid_world = h_mid_local * tier_scale;
    h_mid_world = max(0.0, h_mid_world);

    let up = normalize(p_start);
    let mu = dot(ray_dir, up);
    let curvature_correction = min(1.0 / max(mu, 0.001), 40.0); // Keep the stability clamp

    let optical_length_r = exp(-h_mid_world / h_rayleigh) * path_length * curvature_correction;
    let optical_length_m = exp(-h_mid_world / h_mie) * path_length * curvature_correction;

    return vec2<f32>(optical_length_r, optical_length_m);
}

fn get_sky_color(
    ray: Ray,
    planet_radius_local: f32,
    physical_radius_world: f32,
    planet_center_tier: vec3<f32>,
    light_dir: vec3<f32>, light_color: vec3<f32>,
    scene: SceneUniforms, max_dist_local: f32,
    shadow_casters_in: ptr<storage, array<Sphere>, read>,
    shadow_params_in: ShadowParams,
    ignore_pos_world: vec3<f32>
) -> AtmosphereOutput {
    let params = get_earth_atmosphere();
    let tier_scale = scene.tier_scale_and_pad.x;
    let atmosphere_radius_local = planet_radius_local * ATMOSPHERE_RADIUS_SCALE;

    let r0 = ray.origin;
    let rd = ray.direction;

    let t_atmos = ray_sphere_intersect(r0, rd, atmosphere_radius_local);
    let t_planet = ray_sphere_intersect(r0, rd, planet_radius_local);
    
    let t_start = max(0.0, t_atmos.x);
    var t_end = min(t_atmos.y, max_dist_local);
    if (t_planet.x > t_start) { t_end = min(t_end, t_planet.x); }
    
    let ray_length = t_end - t_start;
    if (ray_length <= 0.0) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0), 0.0); }

    let step_size = ray_length / f32(NUM_IN_SCATTER_SAMPLES);
    var transmittance = vec3<f32>(1.0);
    var scattered_light = vec3<f32>(0.0);

    for (var i = 0u; i < NUM_IN_SCATTER_SAMPLES; i = i + 1u) {
        let p_sample = r0 + rd * (t_start + (f32(i) + 0.5) * step_size);
        if (length(p_sample) < planet_radius_local) { continue; }
        
        let p_world_sq = dot(p_sample, p_sample) * tier_scale * tier_scale;
        let h = (p_world_sq - physical_radius_world*physical_radius_world) / (2.0 * physical_radius_world);

        if (h > 0.0) {
            let density_r = exp(-h / params.h_rayleigh);
            let density_m = exp(-h / params.h_mie);
            let scattering_coeffs = params.beta_rayleigh * density_r + params.beta_mie * density_m;
            var sun_transmittance = vec3<f32>(0.0);
            
            let p_sample_tier = p_sample + planet_center_tier;
            let p_world = p_sample_tier * tier_scale;
            let up_dir = normalize(p_sample);
            let eclipse_shadow_ray = Ray(p_world + up_dir * 0.2, light_dir);
            let eclipse_rec = hit_spheres(eclipse_shadow_ray, createInterval(0.001, INFINITY), shadow_casters_in, shadow_params_in.count, ignore_pos_world);

            if (eclipse_rec.hit && dot(eclipse_rec.emissive, eclipse_rec.emissive) < 0.1) { continue; }

            let optical_lengths = calculate_optical_length_to_space(
                p_sample, light_dir, 
                planet_radius_local,
                atmosphere_radius_local,
                params.h_rayleigh, params.h_mie,
                tier_scale
            );
            if (optical_lengths.x < 1e9) {
                let optical_depth_sun = params.beta_rayleigh * optical_lengths.x + params.beta_mie * optical_lengths.y;
                sun_transmittance = exp(-optical_depth_sun * tier_scale);
            }

            let mu = dot(rd, light_dir);
            let g = params.g_mie;
            let phase_r = 3.0 / (16.0 * PI) * (1.0 + mu*mu);
            let phase_m = 3.0 / (8.0 * PI) * ((1.0-g*g)*(1.0+mu*mu)) / ((2.0+g*g)*pow(1.0+g*g-2.0*g*mu, 1.5));
            let in_scatter = (params.beta_rayleigh * density_r * phase_r) + (params.beta_mie * density_m * phase_m);
				
			scattered_light += transmittance * in_scatter * sun_transmittance * step_size * tier_scale;
			transmittance *= exp(-scattering_coeffs * step_size * tier_scale);
        }
    }
    
    let E_sun = 20.0;
    let final_color = scattered_light * E_sun * light_color;
    
    return AtmosphereOutput(final_color, transmittance, 1.0 - transmittance.g);
}
// Functions for calculating atmospheric in-scattering and transmittance.
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

// Cheap 32-bit integer hash to generate a stable jitter in [0,1)
fn hash11(n: f32) -> f32 {
    let x = fract(n * 0.1031);
    let y = x * (x + 33.33);
    return fract((x + y) * (y + x));
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
    world_scale: f32
) -> vec2<f32> {
    let t_atmos_exit = ray_sphere_intersect(p_start, ray_dir, atmosphere_radius).y;
    if (t_atmos_exit <= 0.0) {
        return vec2<f32>(0.0);
    }
    let t_planet = ray_sphere_intersect(p_start, ray_dir, planet_radius).x;
    if (t_planet > 0.0 && t_planet < t_atmos_exit) {
        return vec2<f32>(1e9);
    }

    let num_samples = NUM_IN_SCATTER_SAMPLES;
    let step_size = t_atmos_exit / f32(num_samples);
    var tau_r = 0.0;
    var tau_m = 0.0;
    for (var i = 0u; i < num_samples; i = i + 1u) {
        let t = (f32(i) + 0.5) * step_size;
        let p = p_start + ray_dir * t;
        let h_world = max(0.0, (length(p) - planet_radius) * world_scale);
        let rho_r = exp(-h_world / h_rayleigh);
        let rho_m = exp(-h_world / h_mie);
        tau_r += rho_r * step_size;
        tau_m += rho_m * step_size;
    }
    return vec2<f32>(tau_r, tau_m);
}

fn get_sky_color(
    ray: Ray,
    planet_radius_local: f32,
    physical_radius_world: f32,
    planet_center_tier: vec3<f32>,
    light_dir: vec3<f32>, light_color: vec3<f32>,
    scene: SceneUniforms, max_dist_local: f32,
    shadow_casters_in: ptr<storage, array<Sphere>, read>,
    shadow_params_in: ShadowParams, camera: CameraUniforms,
    ignore_pos_world: vec3<f32>
) -> AtmosphereOutput {
    let params = get_earth_atmosphere();
    let world_scale = scene.scale_and_flags.x;
    let atmosphere_radius_local = planet_radius_local * ATMOSPHERE_RADIUS_SCALE;

    let r0 = ray.origin;
    let rd = ray.direction;

    let t_atmos = ray_sphere_intersect(r0, rd, atmosphere_radius_local);
    let t_planet = ray_sphere_intersect(r0, rd, planet_radius_local);
    let t_start = max(0.0, t_atmos.x);
    var t_end = min(t_atmos.y, max_dist_local);
    if (t_planet.x > t_start) { t_end = min(t_end, t_planet.x); }

    let num_samples = NUM_IN_SCATTER_SAMPLES;
    let step_size = (t_end - t_start) / f32(num_samples);
    var transmittance = vec3<f32>(1.0);
    var scattered_light = vec3<f32>(0.0);

    for (var i = 0u; i < num_samples; i = i + 1u) {
        let p_sample = r0 + rd * (t_start + (f32(i) + 0.5) * step_size);
        if (length(p_sample) < planet_radius_local) { continue; }
        let h = max(0.0, (length(p_sample) - planet_radius_local) * world_scale);

        if (h > 0.0) {
            let density_r = exp(-h / params.h_rayleigh);
            let density_m = exp(-h / params.h_mie);
            let scattering_coeffs = params.beta_rayleigh * density_r + params.beta_mie * density_m;
            var sun_transmittance = vec3<f32>(0.0);

            let up_dir = normalize(p_sample);
            let eclipse_shadow_ray = Ray(p_sample + up_dir * 0.2, light_dir);
            let eclipse_rec = hit_spheres_shadow(eclipse_shadow_ray, createInterval(0.001, INFINITY), shadow_casters_in, shadow_params_in.count, ignore_pos_world, camera);

            if (eclipse_rec.hit && dot(eclipse_rec.emissive, eclipse_rec.emissive) < 0.1) { continue; }

            let optical_lengths = calculate_optical_length_to_space(
                p_sample, light_dir,
                planet_radius_local, atmosphere_radius_local,
                params.h_rayleigh, params.h_mie,
                world_scale
            );
            let optical_depth_sun = params.beta_rayleigh * optical_lengths.x + params.beta_mie * optical_lengths.y;
            sun_transmittance = exp(-optical_depth_sun * world_scale);

            let mu = dot(rd, light_dir);
            let g = params.g_mie;
            let phase_r = 3.0 / (16.0 * PI) * (1.0 + mu*mu);
            let phase_m = 3.0 / (8.0 * PI) * ((1.0-g*g)*(1.0+mu*mu)) / ((2.0+g*g)*pow(1.0+g*g-2.0*g*mu, 1.5));
            let in_scatter = (params.beta_rayleigh * density_r * phase_r) + (params.beta_mie * density_m * phase_m);
                
            scattered_light += transmittance * in_scatter * sun_transmittance * step_size * world_scale;
            transmittance *= exp(-scattering_coeffs * step_size * world_scale);
        }
    }
    
    let E_sun = 20.0;
    let final_color = scattered_light * E_sun * light_color;
    
    return AtmosphereOutput(final_color, transmittance, 1.0 - transmittance.g);
}



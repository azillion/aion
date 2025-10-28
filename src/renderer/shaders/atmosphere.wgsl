const ATMOSPHERE_RADIUS_SCALE: f32 = 1.025;
const NUM_OUT_SCATTER_SAMPLES: u32 = 8;
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
    // These coefficients are chosen for a clear Earth-like sky.
    // beta_mie is kept low to avoid excessive haze.
    return AtmosphereParams(
        6371.0, 6371.0 * ATMOSPHERE_RADIUS_SCALE,
        // Beta coefficients include 1/wavelength^4 dependency for Rayleigh
        vec3<f32>(5.8e-3, 13.5e-3, 33.1e-3), 8.5,
        vec3<f32>(2.1e-3), 1.2,
        0.76 // Positive for forward scattering
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

fn get_sky_color(
    ray: Ray,
    planet_radius_local: f32,
    physical_radius_world: f32,
    light_dir: vec3<f32>, light_color: vec3<f32>,
    scene: SceneUniforms, max_dist_local: f32
) -> AtmosphereOutput {
    let params = get_earth_atmosphere();
    let tier_scale = scene.tier_scale_and_pad.x;
    let atmosphere_radius_local = planet_radius_local * ATMOSPHERE_RADIUS_SCALE;

    let r0 = ray.origin;
    let rd = ray.direction;

    // --- 1. Calculate ray segment within atmosphere ---
    let t_atmos = ray_sphere_intersect(r0, rd, atmosphere_radius_local);
    let t_planet = ray_sphere_intersect(r0, rd, planet_radius_local);
    
    let t_start = max(0.0, t_atmos.x);
    var t_end = min(t_atmos.y, max_dist_local);
    if (t_planet.x > t_start) { t_end = min(t_end, t_planet.x); }
    
    let ray_length = t_end - t_start;
    if (ray_length <= 0.0) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0), 0.0); }

    // --- 2. Sampling ---
    // Fixed sample counts for stability and quality.
    let num_in_scatter_samples = 16;
    let num_out_scatter_samples = 8;

    // --- 3. Initialize accumulators ---
    let step_size = ray_length / f32(num_in_scatter_samples);
    var transmittance = vec3<f32>(1.0);
    var scattered_light = vec3<f32>(0.0);

    // --- 4. March along the view ray ---
    for (var i = 0; i < num_in_scatter_samples; i = i + 1) {
        let p_sample = r0 + rd * (t_start + (f32(i) + 0.5) * step_size);
        if (length(p_sample) < planet_radius_local) { continue; }
        // Numerically stable altitude: (|p|^2 - R^2) / (2R), with |p| in world units.
        let p_world_sq = dot(p_sample, p_sample) * tier_scale * tier_scale;
        let h = (p_world_sq - physical_radius_world*physical_radius_world) / (2.0 * physical_radius_world);

        if (h > 0.0) {
            let density_r = exp(-h / params.h_rayleigh);
            let density_m = exp(-h / params.h_mie);
            let scattering_coeffs = params.beta_rayleigh * density_r + params.beta_mie * density_m;

            // --- 4a. Sun visibility and out-scattering along the sun ray ---
            // Compute transmittance only if the sun ray does NOT intersect the planet.
            var sun_transmittance = vec3<f32>(0.0);
            let t_to_planet_surface = ray_sphere_intersect(p_sample, light_dir, planet_radius_local).x;
            if (t_to_planet_surface < 0.0 || t_to_planet_surface > 1.0e4) {
                let t_sun_atmos = ray_sphere_intersect(p_sample, light_dir, atmosphere_radius_local).y;
                var optical_depth_sun = vec3<f32>(0.0);
					if (t_sun_atmos > 0.0) {
						let sun_step_size = t_sun_atmos / f32(num_out_scatter_samples);
						for (var j = 0; j < num_out_scatter_samples; j = j + 1) {
                        let p_sun_sample = p_sample + light_dir * (f32(j) + 0.5) * sun_step_size;
                        // Stable altitude for sun-ray samples as well
                        let p_sun_world_sq = dot(p_sun_sample, p_sun_sample) * tier_scale * tier_scale;
                        let h_sun = (p_sun_world_sq - physical_radius_world*physical_radius_world) / (2.0 * physical_radius_world);
							if (h_sun > 0.0) {
								optical_depth_sun += (
									params.beta_rayleigh * exp(-h_sun / params.h_rayleigh) +
									params.beta_mie * exp(-h_sun / params.h_mie)
								) * sun_step_size * tier_scale;
							}
						}
					}
                sun_transmittance = exp(-optical_depth_sun);
            }

            // --- 4b. Accumulate in-scattered light ---
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

// --- Dummy functions for compatibility with tierCompute.wgsl ---
// These ensure the shader compiles but are not used by the main logic.
struct AtmosphereParamsDummy { dummy: f32 };
fn get_optical_length(p_start: vec3<f32>, ray_dir: vec3<f32>, length: f32, params_dummy: AtmosphereParams, tier_scale: f32) -> vec2<f32> {
    return vec2<f32>(0.0);
}
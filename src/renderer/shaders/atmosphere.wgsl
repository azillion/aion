const PI: f32 = 3.1415926535;

// --- Constants based on the GPU Gems 2 chapter ---
// All calculations are done in a "unit space" where the planet's radius is 1.0.
const R_INNER: f32 = 1.0;
const R_ATMOS: f32 = 1.025; // Atmosphere is ~2.5% of the planet's radius
const ATMOS_THICKNESS: f32 = R_ATMOS - R_INNER; // Normalization denominator for heights

const NUM_SAMPLES: i32 = 16; // Number of samples along the view ray

// Physical/scattering constants
const K_R: f32 = 0.0025; // Rayleigh scattering coefficient
const K_M: f32 = 0.0010; // Mie scattering coefficient
const E_SUN: f32 = 20.0; // Sun intensity (artistic)

// Visible wavelengths in micrometers (approx.) for RGB
const WAVELENGTHS: vec3<f32> = vec3<f32>(0.650, 0.570, 0.475);
const INV_WAVELENGTH4: vec3<f32> = vec3<f32>(
    pow(WAVELENGTHS.x, -4.0),
    pow(WAVELENGTHS.y, -4.0),
    pow(WAVELENGTHS.z, -4.0),
);

// Mie phase asymmetry factor (high magnitude → strong forward scattering)
const MIE_G: f32 = -0.990;

// Scale height (H0 in the article). The atmosphere's density is exp(-h/H0).
// The article commonly uses ~0.25 for nice falloff in unit space.
const SCALE_DEPTH: f32 = 0.25;
const INV_SCALE_DEPTH: f32 = 1.0 / SCALE_DEPTH;

struct AtmosphereOutput {
	color: vec3<f32>,
	transmittance: vec3<f32>,
};

// --- Helper Functions ---

// Ray-sphere intersection in unit space; returns near/far t.
fn ray_vs_sphere(p: vec3<f32>, dir: vec3<f32>, r: f32) -> vec2<f32> {
	let b = dot(p, dir);
	let c = dot(p, p) - r * r;
	var d = b * b - c;
	if (d < 0.0) { return vec2<f32>(1e30, -1e30); }
	d = sqrt(d);
	return vec2<f32>(-b - d, -b + d);
}

// The "scale" function from GPU Gems. Approximates optical depth based on angle.
// `mu` is the cosine of the angle between the ray and the local vertical (up).
fn scale(mu: f32) -> f32 {
	let x = 1.0 - mu;
	// This polynomial is a curve-fit approximation of the optical depth integral based on angle.
	// The `exp()` wrapper found in the original code is incorrect and leads to astronomical values,
	// causing floating point overflow, NaN results, and black pixels.
	// The result is also clamped to prevent negative optical depth.
	let polynomial = -0.00287 + x*(0.459 + x*(3.83 + x*(-6.80 + x*5.25)));
	return max(0.0, polynomial);
}

// Phase functions
fn phase_ray(mu: f32) -> f32 {
	return (3.0 / (16.0 * PI)) * (1.0 + mu * mu);
}

fn phase_mie(mu: f32) -> f32 {
	let g = MIE_G;
	let gg = g * g;
	let a = (1.0 - gg) * (1.0 + mu*mu);
	let b = (2.0 + gg) * pow(1.0 + gg - 2.0 * g * mu, 1.5);
	return (3.0 / (8.0 * PI)) * (a / b);
}

// --- Dedicated sunlight transmittance from surface to sun ---
// Computes how the atmosphere filters sunlight traveling from a surface point toward the sun.
// Returns 0 when the sun is below the local horizon (blocked by the planet).
fn get_sunlight_transmittance(
	surface_pos_unit: vec3<f32>,
	light_dir: vec3<f32>
) -> vec3<f32> {
	let up = normalize(surface_pos_unit);
	let ld = normalize(light_dir);
	// Stable horizon check: sun is below horizon if it's pointing below the local up
	let sun_is_below_horizon = dot(up, ld) < 0.0;
	if (sun_is_below_horizon) { return vec3<f32>(0.0); }

	let surface_height = length(surface_pos_unit);
	let h_norm = (surface_height - R_INNER) / ATMOS_THICKNESS;
	let depth = exp(-h_norm * INV_SCALE_DEPTH);
	let mu = dot(ld, up);
	let scatter = depth * scale(mu);

	let extinction_ray = K_R * 4.0 * PI * INV_WAVELENGTH4;
	let extinction_mie = K_M * 4.0 * PI;
	let extinction = extinction_ray + vec3<f32>(extinction_mie);

	return exp(-scatter * extinction);
}

// --- Dedicated fast transmittance function ---
// Computes only the transmittance along a ray in unit space using the scale function.
fn get_transmittance(
	ray_origin_unit: vec3<f32>,
	ray_dir: vec3<f32>
) -> vec3<f32> {
	let cam_height = length(ray_origin_unit);

	let e = ray_vs_sphere(ray_origin_unit, ray_dir, R_ATMOS);
	if (e.x > e.y) { return vec3<f32>(1.0); }

	let h_norm = (cam_height - R_INNER) / ATMOS_THICKNESS;
	let depth = exp(-h_norm * INV_SCALE_DEPTH);
	let mu = dot(ray_dir, ray_origin_unit) / cam_height;
	let scatter = depth * scale(mu);

	let extinction_ray = K_R * 4.0 * PI * INV_WAVELENGTH4;
	let extinction_mie = K_M * 4.0 * PI;
	let extinction = extinction_ray + vec3<f32>(extinction_mie);

	return exp(-scatter * extinction);
}

// --- Main API ---
// Calculates in-scattered sky radiance and view-ray transmittance.

fn get_sky_color(
	ray: Ray,
	planet_center_tier: vec3<f32>,
	planet_radius_tier: f32,
	light_dir: vec3<f32>,
	scene: SceneUniforms
) -> AtmosphereOutput {
	// 1) Convert problem into unit space where the planet radius is 1.0
	let o = (ray.origin - planet_center_tier) / planet_radius_tier;
	let dir = ray.direction;
	let cam_height = length(o);

	// 2) Intersect view ray with atmospheric shell, clip against solid planet
	var e = ray_vs_sphere(o, dir, R_ATMOS);
	if (e.x > e.y) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0)); }
	let inner = ray_vs_sphere(o, dir, R_INNER);
	e.y = min(e.y, inner.x);
	if (e.x > e.y) { return AtmosphereOutput(vec3<f32>(0.0), vec3<f32>(1.0)); }
	e.x = max(e.x, 0.0);

	// 3) Precompute terms used in the sampling loop (use normalized height)
	let h_norm_start = (cam_height - R_INNER) / ATMOS_THICKNESS;
	let start_depth = exp(-h_norm_start * INV_SCALE_DEPTH);
	let start_angle_mu = dot(o, dir) / cam_height;
	let start_offset = start_depth * scale(start_angle_mu);

	let segment_len = (e.y - e.x) / f32(NUM_SAMPLES);
	var sample_pos = o + dir * (e.x + segment_len * 0.5);

	var front_color = vec3<f32>(0.0);

	let extinction = INV_WAVELENGTH4 * (K_R * 4.0 * PI) + vec3<f32>(K_M * 4.0 * PI);
	let light_dir_n = normalize(light_dir);

	// 4) Single outer loop along the view ray, inner depth via scale function
	for (var i = 0; i < NUM_SAMPLES; i = i + 1) {
		let sample_height = length(sample_pos);
		if (sample_height < R_INNER) { sample_pos += dir * segment_len; continue; }

		// Density must use normalized height in [0, 1]
		let h_norm_sample = (sample_height - R_INNER) / ATMOS_THICKNESS;
		let local_density = exp(-h_norm_sample * INV_SCALE_DEPTH);
		let light_angle_mu = dot(light_dir_n, sample_pos) / sample_height;
		let camera_angle_mu = dot(dir, sample_pos) / sample_height;
		let scatter = max(0.0, start_offset + local_density * (scale(light_angle_mu) - scale(camera_angle_mu)));

		let att = exp(-scatter * extinction);
		front_color += att * (local_density * segment_len);
		sample_pos += dir * segment_len;
	}

	// 5) Final color using phase functions
	let mu = dot(dir, light_dir_n);
	let rayleigh_color = front_color * (INV_WAVELENGTH4 * (K_R * E_SUN)) * phase_ray(mu);
	let mie_color = front_color * vec3<f32>(K_M * E_SUN) * phase_mie(mu);
	let final_color = (rayleigh_color + mie_color) * scene.dominant_light_color_and_debug.xyz;

	// 6) View ray transmittance for surface lighting (use normalized height)
	let h_norm_view = (cam_height - R_INNER) / ATMOS_THICKNESS;
	let view_depth = exp(-h_norm_view * INV_SCALE_DEPTH);
	let view_angle_mu = dot(o, dir) / cam_height;
	let view_scatter = view_depth * scale(view_angle_mu);
	let view_transmittance = exp(-view_scatter * extinction);

	return AtmosphereOutput(final_color, view_transmittance);
}
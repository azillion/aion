// Functions for procedural terrain generation and ray-marching Signed Distance Fields (SDFs).

// Domain warping function for more interesting shapes.
fn warp(p: vec3<f32>, seed: f32) -> vec3<f32> {
    let q2 = vec2<f32>(
        snoise(p + vec3<f32>(0.0, 0.0, seed)),
        snoise(p + vec3<f32>(5.2, 1.3, seed))
    );
    return p + 0.4 * vec3<f32>(q2.x, q2.y, q2.x);
}

// Fractal Brownian Motion (FBM) using simplex noise
fn fbm(p: vec3<f32>, octaves: u32) -> f32 {
    var total = 0.0;
    var frequency = 1.0;
    var amplitude = 0.5;
    let lacunarity = 2.0;
    let persistence = 0.5;

    for (var i = 0u; i < octaves; i = i + 1u) {
        total += snoise(p * frequency) * amplitude;
        frequency *= lacunarity;
        amplitude *= persistence;
    }
    return total;
}

// Single-pass fractional FBM to reduce snoise calls
fn fbm_frac(p: vec3<f32>, octaves_f: f32) -> f32 {
    let o_floor = floor(octaves_f);
    let oi = u32(o_floor);
    let frac = octaves_f - o_floor;
    var total = 0.0;
    var freq  = 1.0;
    var amp   = 0.5;
    for (var i = 0u; i < oi; i = i + 1u) {
        total += snoise(p * freq) * amp;
        freq  *= 2.0;
        amp   *= 0.5;
    }
    if (frac > 0.0) {
        total += snoise(p * freq) * amp * frac;
    }
    return total;
}

fn remap01_any(a: f32, b: f32, x: f32) -> f32 {
	let t = clamp((x - a) / (b - a), 0.0, 1.0);
	return select(t, 1.0 - t, a > b);
}

fn remap01(a: f32, b: f32, x: f32) -> f32 { return clamp((x - a) / (b - a), 0.0, 1.0); }
fn smin(a: f32, b: f32, k: f32) -> f32 {
    let h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return mix(b, a, h) - k * h * (1.0 - h);
}

// Calculates signed terrain height (km) using FBM with stable input domain
fn h_noise(p_local: vec3<f32>, params: TerrainUniforms, dist_to_surface: f32, camera: CameraUniforms, scene: SceneUniforms) -> f32 {
    // Stable domain: normalize local position to get a precise direction
    let dir = normalize(p_local);

    // Smooth LOD: vary octaves with distance to surface
    let min_oct = 3.0;
    let max_oct = 7.0;
    let lod_start = params.base_radius * 2.0;
    let lod_end = params.base_radius * 0.1;
    let lod_t = remap01_any(lod_start, lod_end, dist_to_surface);
    let octaves_f = mix(min_oct, max_oct, lod_t);

    // Domain warp on stable direction
    let base_frequency = 2.0;
    let p_warped = warp(dir * base_frequency, params.seed);

    // FBM in [-1,1]
    let n = fbm_frac(p_warped, octaves_f);

    // Slope-aware amplitude scaling to self-guard when octaves increase
    let s_max = 0.8;
    let lipschitz_fbm = base_frequency * octaves_f;
    let amp_scale = s_max / max(1.0, lipschitz_fbm);

    // Map to [0,1] and scale to kilometers
    let height_range = params.max_height * params.base_radius * amp_scale;
    return (n * 0.5 + 0.5) * height_range;
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
    let tier_scale = scene.scale_and_flags.x;

    var up_vec = vec3<f32>(0.0, 1.0, 0.0);
    if (abs(dot(analytic_normal, up_vec)) > 0.999) { up_vec = vec3<f32>(1.0, 0.0, 0.0); }
    let tangent = normalize(cross(analytic_normal, up_vec));
    let bitangent = normalize(cross(analytic_normal, tangent));

    let h0 = h_noise(p_on_sphere - planet_center, params, dist_to_surface, camera, scene);
    let h1 = h_noise(p_on_sphere + tangent * eps - planet_center, params, dist_to_surface + eps, camera, scene);
    let h2 = h_noise(p_on_sphere + bitangent * eps - planet_center, params, dist_to_surface + eps, camera, scene);

    let final_normal = normalize(
        analytic_normal - (tangent * (h1 - h0) / eps + bitangent * (h2 - h0) / eps) * (1.0 / tier_scale)
    );
    return final_normal;
}

fn get_sdf_normal(
    p_local: vec3<f32>,
    params: TerrainUniforms,
    t: f32,
    camera: CameraUniforms,
    scene: SceneUniforms,
    lod_dist: f32
) -> vec3<f32> {
    // Stable central differencing with epsilon tied to marched distance
    let eps = max(0.001, t * 0.0001);
    let dx = vec3<f32>(eps, 0.0, 0.0);
    let dy = vec3<f32>(0.0, eps, 0.0);
    let dz = vec3<f32>(0.0, 0.0, eps);

    return normalize(vec3<f32>(
        dWorld(p_local + dx, params, t, camera, scene, lod_dist) - dWorld(p_local - dx, params, t, camera, scene, lod_dist),
        dWorld(p_local + dy, params, t, camera, scene, lod_dist) - dWorld(p_local - dy, params, t, camera, scene, lod_dist),
        dWorld(p_local + dz, params, t, camera, scene, lod_dist) - dWorld(p_local - dz, params, t, camera, scene, lod_dist)
    ));
}

fn dWorld(
    p_local: vec3<f32>,
    params: TerrainUniforms,
    dist_marched: f32,
    camera: CameraUniforms,
    scene: SceneUniforms,
    lod_dist: f32
) -> f32 {
    let tier_scale = scene.scale_and_flags.x;
    let inv_tier = 1.0 / tier_scale;
    let len_p = length(p_local);
    let h = h_noise(p_local, params, lod_dist, camera, scene);
    let h_scaled = h * inv_tier;
    let R_scaled = params.base_radius * inv_tier;
    let d_terrain = len_p - (R_scaled + h_scaled);
    
    // Dynamic water level from simulation grid
    let dir = normalize(p_local);
    let water_uv_face = directionToCubeUV(dir);
    let dims = textureDimensions(water_state, 0);
    let water_sample = textureLoad(water_state, vec2<i32>(water_uv_face.xy * vec2<f32>(dims)), i32(water_uv_face.z), 0);
    let water_depth = water_sample.r; // km
    let water_radius_scaled = (params.base_radius + params.sea_level + water_depth) * inv_tier;
    let d_ocean = len_p - water_radius_scaled;
    
    let k = 0.01 * tier_scale;
    let m = min(d_terrain, d_ocean);
    let sd = smin(d_terrain, d_ocean, k);
    return select(m, sd, abs(d_terrain - d_ocean) < 2.0 * k);
}

// Local helper: ray-sphere intersection in planet-local space
fn ray_sphere_intersect_local(r0: vec3<f32>, rd: vec3<f32>, sr: f32) -> vec2<f32> {
    let b = dot(rd, r0);
    let c = dot(r0, r0) - (sr * sr);
    let d = b * b - c;
    if (d < 0.0) { return vec2<f32>(1e5, -1e5); }
    return vec2<f32>(-b - sqrt(d), -b + sqrt(d));
}

// High-precision ray type using f64 emulation (high/low pairs per component)
struct Ray_f64 { origin_h: vec3<f32>, origin_l: vec3<f32>, direction: vec3<f32> };

fn rayAt_f64(ray: Ray_f64, t: f32) -> vec3<f32> {
    let tdir = ray.direction * t;
    let px = add_f64_f32(vec2<f32>(ray.origin_h.x, ray.origin_l.x), tdir.x);
    let py = add_f64_f32(vec2<f32>(ray.origin_h.y, ray.origin_l.y), tdir.y);
    let pz = add_f64_f32(vec2<f32>(ray.origin_h.z, ray.origin_l.z), tdir.z);
    return vec3<f32>(px.x + px.y, py.x + py.y, pz.x + pz.y);
}

fn ray_march(
	ray_local: Ray_f64,
    max_dist: f32,
    params: TerrainUniforms,
    camera: CameraUniforms,
    scene: SceneUniforms,
    lod_dist: f32
) -> HitRecord {
	var rec: HitRecord; rec.hit = false;

	// Compute a stable local starting point by intersecting a conservative bounding shell
	let tier_scale = scene.scale_and_flags.x;
	let R_scaled = params.base_radius / tier_scale;
	let max_height_scaled = (params.max_height * params.base_radius) / tier_scale;
	let water_radius_scaled = (params.base_radius + params.sea_level) / tier_scale;
	let shell_radius = max(R_scaled + max_height_scaled * 2.0, water_radius_scaled);

	let origin_f32 = ray_local.origin_h + ray_local.origin_l;
	let t_shell = ray_sphere_intersect_local(origin_f32, ray_local.direction, shell_radius);
	let t_start = max(0.0, t_shell.x);
	let t_end = min(max_dist, t_shell.y);
	if (t_end <= t_start) { return rec; }

	var t = t_start;
	for (var i = 0; i < 64; i = i + 1) {
		let p_local = rayAt_f64(ray_local, t);
		let d_scaled = dWorld(p_local, params, t, camera, scene, lod_dist);
		let d = d_scaled / tier_scale;

		// Relative threshold tied to camera distance to avoid misses at limb
		let threshold = 0.0001 * t;
		if (d < threshold) {
			rec.hit = true;
			rec.t = t;
			rec.p = p_local;
			let outward_normal = get_sdf_normal(p_local, params, t, camera, scene, lod_dist);
			rec.front_face = dot(ray_local.direction, outward_normal) < 0.0;
			rec.normal = outward_normal;
			return rec;
		}

		t = t + d;
		if (t > t_end) { break; }
	}

	return rec;
}


// (Material and shading helpers moved to shading.wgsl)

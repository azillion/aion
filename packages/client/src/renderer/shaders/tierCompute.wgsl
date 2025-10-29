#include "camera.wgsl"
#include "sceneUniforms.wgsl"
#include "noise.wgsl"


const INFINITY: f32 = 1e38;
const PI: f32 = 3.1415926535;

fn degreesToRadians(degrees: f32) -> f32 { return degrees * PI / 180.0; }

// Converts a normalized 3D position on a sphere to 2D equirectangular UV coordinates.
fn worldPosToUV(p_norm: vec3<f32>) -> vec2<f32> {
    let phi = atan2(p_norm.y, p_norm.x);
    let theta = acos(clamp(p_norm.z, -1.0, 1.0));
    let u = (phi + PI) / (2.0 * PI);
    let v = theta / PI;
    return vec2<f32>(u, v);
}

struct Interval { minI: f32, maxI: f32 };
fn createInterval(min: f32, max: f32) -> Interval { return Interval(min, max); }
fn intervalSurrounds(i: Interval, x: f32) -> bool { return i.minI < x && x < i.maxI; }

struct Sphere {
	pos_and_radius: vec4<f32>,              // .xyz = position, .w = geometric radius
	albedo_and_atmos_flag: vec4<f32>,       // .xyz = albedo, .w = has_atmosphere_flag
	emissive_and_terrain_flag: vec4<f32>,   // .xyz = emissive, .w = has_terrain_flag
	ref_idx_opacity_pad: vec4<f32>,         // .x = ref_idx, .y = opacity, .zw unused (placeholder for alignment)
	terrain_params: vec4<f32>,              // .x = base_radius, .y = sea_level, .z = max_height, .w = seed
	padding_vec: vec4<f32>,                 // Padding to maintain 24-float stride
}

struct HitRecord {
    p: vec3<f32>, normal: vec3<f32>, t: f32, hit: bool, front_face: bool,
    albedo: vec3<f32>, emissive: vec3<f32>, object_index: u32,
}

struct Ray { origin: vec3<f32>, direction: vec3<f32> }

@group(0) @binding(6) var<storage, read> shadow_casters: array<Sphere>;
struct ShadowParams { count: u32 };
@group(0) @binding(7) var<uniform> shadowParams: ShadowParams;

fn rayAt(ray: Ray, t: f32) -> vec3<f32> { return ray.origin + ray.direction * t; }

fn hit_sphere(sphere: Sphere, oc: vec3<f32>, r: Ray, ray_t: Interval, index: u32) -> HitRecord {
    var rec: HitRecord; rec.hit = false;
    let radius = sphere.pos_and_radius.w;
    let a = dot(r.direction, r.direction); let h = dot(r.direction, oc); let c = dot(oc, oc) - radius * radius;
    let discriminant = h * h - a * c; if (discriminant < 0.0) { return rec; }
    let sqrtd = sqrt(discriminant);
    var root = (h - sqrtd) / a; if (!intervalSurrounds(ray_t, root)) { root = (h + sqrtd) / a; if (!intervalSurrounds(ray_t, root)) { return rec; } }
    rec.t = root; rec.p = rayAt(r, rec.t);
    let outward_normal = (rec.p - sphere.pos_and_radius.xyz) / radius;
    rec.front_face = dot(r.direction, outward_normal) < 0.0;
    rec.normal = outward_normal;
    rec.albedo = sphere.albedo_and_atmos_flag.xyz; rec.emissive = sphere.emissive_and_terrain_flag.xyz; rec.object_index = index; rec.hit = true; return rec;
}

// Performs a numerically stable check against a list of spheres and returns the closest hit.
fn hit_spheres(r: Ray, ray_t: Interval, spheres_in: ptr<storage, array<Sphere>, read>, sphere_count: u32, ignore_pos: vec3<f32>) -> HitRecord {
    var closest_so_far = ray_t.maxI;
    var rec: HitRecord;
    rec.hit = false;
    for (var i = 0u; i < sphere_count; i++) {
        let sphere = (*spheres_in)[i];
        // Ignore if the sphere center is very close to the position we want to ignore (e.g., self-shadowing).
        if (distance(sphere.pos_and_radius.xyz, ignore_pos) < 0.1) { continue; }

        // CRITICAL FIX: The subtraction of two large world-space vectors (sphere center and ray origin)
        // is numerically unstable with f32. By performing it here once and passing the result (`oc`)
        // to the intersection function, we perform the sensitive operation in the calling context
        // and keep the core intersection logic stable and reusable.
        let oc = sphere.pos_and_radius.xyz - r.origin;
        let sphere_rec = hit_sphere(sphere, oc, r, Interval(ray_t.minI, closest_so_far), i);
        if (sphere_rec.hit) {
            closest_so_far = sphere_rec.t;
            rec = sphere_rec;
        }
    }
    return rec;
}

#include "atmosphere.wgsl"
 

struct Camera { origin: vec3<f32>, lower_left_corner: vec3<f32>, horizontal: vec3<f32>, vertical: vec3<f32> }

fn createCamera(aspect_ratio: f32) -> Camera {
    let vfov = 25.0; let lookfrom = vec3<f32>(0.0, 0.0, 0.0); let vup = camera.up;
    let theta = degreesToRadians(vfov); let h = tan(theta / 2.0); let viewport_height = 2.0 * h; let viewport_width = aspect_ratio * viewport_height;
    let w = -camera.forward; let u = normalize(cross(vup, w)); let v = cross(w, u);
    let origin = lookfrom; let horizontal = viewport_width * u; let vertical = viewport_height * v; let lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - w;
    return Camera(origin, lower_left_corner, horizontal, vertical);
}

fn getRay(camera: Camera, s: f32, t: f32) -> Ray {
    return Ray(camera.origin, normalize(camera.lower_left_corner + s*camera.horizontal + t*camera.vertical - camera.origin));
}

fn shade_planet_surface(rec: HitRecord, sphere: Sphere, ray: Ray, scene: SceneUniforms, camera: CameraUniforms) -> vec3<f32> {
    // let terrain_uniforms = TerrainUniforms( sphere.terrain_params.x, sphere.terrain_params.y, sphere.terrain_params.z, sphere.terrain_params.w );
    // let tier_scale = scene.tier_scale_and_pad.x;
    // let R_world = terrain_uniforms.base_radius;
	// // Use a consistent LOD distance per-planet for this frame to avoid discontinuities
	// let dist_to_center = length(sphere.pos_and_radius.xyz);
	// let lod_dist_world = dist_to_center * tier_scale;
    // let p_local_hit = rec.p - sphere.pos_and_radius.xyz;
	// let h_world = h_noise(normalize(p_local_hit), terrain_uniforms, lod_dist_world, R_world, camera, scene);
    // let is_ocean = h_world < terrain_uniforms.sea_level;

    // var surface_albedo: vec3<f32>; var surface_normal: vec3<f32>;
    // if (is_ocean) {
    //     surface_albedo = get_ocean_material(p_local_hit).albedo;
    //     surface_normal = get_ocean_wave_normal(p_local_hit, normalize(p_local_hit));
    // } else {
	// 	surface_normal = get_terrain_normal_from_heightfield( rec.p, sphere.pos_and_radius.xyz, rec.normal, terrain_uniforms, lod_dist_world, R_world, camera, scene );
    //     if (terrain_uniforms.seed == 1337.0) {
    //         surface_albedo = get_lunar_material(p_local_hit, surface_normal).albedo;
    //     } else {
    //         surface_albedo = get_material(p_local_hit, surface_normal, terrain_uniforms, h_world).albedo;
    //     }
    // }

    // --- DEBUG: Simplified shading for a smooth sphere to isolate atmosphere issues ---
    let surface_albedo = sphere.albedo_and_atmos_flag.xyz;
    let surface_normal = rec.normal;
    let is_ocean = sphere.terrain_params.y > 0.0; // Simple ocean check based on sea level

    let light_dir_to_source = -normalize(scene.dominant_light_direction.xyz);
    let light_color = scene.dominant_light_color_and_debug.xyz;
    let view_dir = normalize(ray.origin - rec.p);
    let tier_scale = scene.tier_scale_and_pad.x;
    let p_world = rec.p * tier_scale;
	// The normal in world space is unchanged under uniform scale; using rec.normal directly is correct.
	let self_pos_world = sphere.pos_and_radius.xyz * tier_scale;
	let shadow_ray = Ray(p_world + surface_normal * 0.2, light_dir_to_source);
	let inter_body_shadow_rec = hit_spheres(shadow_ray, createInterval(0.001, INFINITY), &shadow_casters, shadowParams.count, self_pos_world);
	let is_eclipsed = inter_body_shadow_rec.hit && dot(inter_body_shadow_rec.emissive, inter_body_shadow_rec.emissive) < 0.1;

	// Sunlight transmittance through atmosphere to the surface point
	var sun_color = light_color;
	let has_atmosphere = sphere.albedo_and_atmos_flag.w > 0.5;
	var lit_color = surface_albedo * select(0.0, 0.02, has_atmosphere); // Start with a minimal non-atmospheric ambient

	if (has_atmosphere) {
		// --- UNIFIED LIGHTING MODEL ---
		// The surface is ALWAYS lit by the sky. This is our base ambient term.
		let zenith_dir = surface_normal;
		let sky_ambient = get_sky_color(
			Ray(rec.p - sphere.pos_and_radius.xyz, zenith_dir),
			sphere.pos_and_radius.w, sphere.terrain_params.x,
			sphere.pos_and_radius.xyz, light_dir_to_source, light_color, scene, INFINITY,
			&shadow_casters, shadowParams, self_pos_world
		);
		lit_color = surface_albedo * sky_ambient.in_scattering;

		// The surface is ADDITIONALLY lit by the sun if it's above the horizon.
		let dot_nl = dot(surface_normal, light_dir_to_source);
		if (dot_nl > 0.0 && !is_eclipsed) {
			// Calculate transmittance for the direct sunbeam
			 let params = get_earth_atmosphere();
			 let tier_scale = scene.tier_scale_and_pad.x;
			 let p_start_local = rec.p - sphere.pos_and_radius.xyz;
			 let optical_lengths = calculate_optical_length_to_space(
			 	 p_start_local, light_dir_to_source,
			 	 sphere.pos_and_radius.w, sphere.pos_and_radius.w * ATMOSPHERE_RADIUS_SCALE,
			 	 params.h_rayleigh, params.h_mie, tier_scale
			 );
			 let optical_depth = (params.beta_rayleigh * optical_lengths.x + params.beta_mie * optical_lengths.y) * tier_scale;
			 let transmittance = exp(-optical_depth);
			 let final_sun_color = light_color * transmittance;

			// Add diffuse and specular from the direct sunbeam
			 let diffuse_term = surface_albedo * dot_nl * final_sun_color;
			 let half_vec = normalize(light_dir_to_source + view_dir);
			 let specular_term = pow(max(0.0, dot(surface_normal, half_vec)), 64.0) * select(0.1, 1.0, is_ocean) * final_sun_color;
			 lit_color += (diffuse_term + specular_term);
		}
	}
	
	return lit_color;
}

 

struct ComputeParams { bodyCount: u32 };
@group(0) @binding(0) var<uniform> params: ComputeParams;
@group(0) @binding(1) var<storage, read> spheres: array<Sphere>;
@group(0) @binding(2) var output: texture_storage_2d<rgba16float, write>;
@group(0) @binding(3) var<uniform> camera: CameraUniforms;
@group(0) @binding(4) var depthOut: texture_storage_2d<r32float, write>;
@group(0) @binding(5) var<uniform> scene: SceneUniforms;
@group(0) @binding(8) var gol_texture: texture_2d<u32>;


@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
	let dims = textureDimensions(output);
	if (id.x >= dims.x || id.y >= dims.y) { return; }
	let coords_i = vec2<i32>(id.xy);

	let aspect_ratio = f32(dims.x) / f32(dims.y);
	let cam = createCamera(aspect_ratio);
	let u = (f32(coords_i.x) + 0.5) / f32(dims.x);
	let v = 1.0 - (f32(coords_i.y) + 0.5) / f32(dims.y);
	let ray = getRay(cam, u, v);

	let atmosphereEnabled = scene.tier_scale_and_pad.y > 0.5;

	// --- NEW, HIGH-PERFORMANCE RENDER LOOP ---

	// 1. BROAD PHASE: Find the closest analytical sphere hit. This is very fast.
    var rec = hit_spheres(ray, createInterval(0.001, INFINITY), &spheres, params.bodyCount, vec3<f32>(1e30));

	var final_pixel_color = vec3<f32>(0.0);
	var final_alpha = 0.0;

	let max_dist = select(INFINITY, rec.t, rec.hit);

	var total_transmittance = vec3<f32>(1.0);
	var total_in_scattering = vec3<f32>(0.0);

	if (atmosphereEnabled) {
		for (var i = 0u; i < params.bodyCount; i = i + 1u) {
			if (rec.hit && rec.object_index == i) { continue; }
			let atmos_sphere = spheres[i];
			if (atmos_sphere.albedo_and_atmos_flag.w > 0.5) {
				let self_pos_world = atmos_sphere.pos_and_radius.xyz * scene.tier_scale_and_pad.x;
				let atmos = get_sky_color(
					Ray(ray.origin - atmos_sphere.pos_and_radius.xyz, ray.direction),
					atmos_sphere.pos_and_radius.w, atmos_sphere.terrain_params.x,
					atmos_sphere.pos_and_radius.xyz,
					-normalize(scene.dominant_light_direction.xyz),
					scene.dominant_light_color_and_debug.xyz, //
					scene, max_dist, &shadow_casters, shadowParams, self_pos_world
				);
				total_in_scattering += atmos.in_scattering * total_transmittance;
				total_transmittance *= atmos.transmittance;
				final_alpha = max(final_alpha, atmos.alpha);
			}
		}
	}

	if (rec.hit) {
		let hit_sphere = spheres[rec.object_index];
		var surface_color = vec3<f32>(0.0);

		if (dot(hit_sphere.emissive_and_terrain_flag.xyz, hit_sphere.emissive_and_terrain_flag.xyz) > 0.1) {
			surface_color = hit_sphere.emissive_and_terrain_flag.xyz;
		} else if (atmosphereEnabled && hit_sphere.albedo_and_atmos_flag.w > 0.5) {
			// 1. Get the color of the surface, lit by sunlight that has passed through the atmosphere.
			let ground_color = shade_planet_surface(rec, hit_sphere, ray, scene, camera);

			// 2. Calculate atmospheric scattering along the view ray (from camera to the surface).
			let self_pos_world = hit_sphere.pos_and_radius.xyz * scene.tier_scale_and_pad.x;
			let view_atmos = get_sky_color(
				Ray(ray.origin - hit_sphere.pos_and_radius.xyz, ray.direction),
				hit_sphere.pos_and_radius.w, hit_sphere.terrain_params.x,
				hit_sphere.pos_and_radius.xyz,
				-normalize(scene.dominant_light_direction.xyz),
				scene.dominant_light_color_and_debug.xyz,
				scene, rec.t, &shadow_casters, shadowParams, self_pos_world
			);

			// 3. Combine them: Surface color is attenuated by the view ray, then in-scattered light is added.
			surface_color = ground_color * view_atmos.transmittance + view_atmos.in_scattering;
		} else { // Moon or other simple sphere
			let light_dir_to_source = -normalize(scene.dominant_light_direction.xyz);
			let light_color = scene.dominant_light_color_and_debug.xyz;
			let diffuse_intensity = max(dot(rec.normal, light_dir_to_source), 0.0);
			let p_world = rec.p * scene.tier_scale_and_pad.x;
			let shadow_ray = Ray(p_world + rec.normal * 0.2, light_dir_to_source);
			let self_pos_world_unscaled = hit_sphere.pos_and_radius.xyz * scene.tier_scale_and_pad.x;
			let shadow_rec = hit_spheres(shadow_ray, createInterval(0.001, INFINITY), &shadow_casters, shadowParams.count, self_pos_world_unscaled);
			var shadow_multiplier = 1.0;
            if (shadow_rec.hit && dot(shadow_rec.emissive, shadow_rec.emissive) < 0.1) { shadow_multiplier = 0.01; }
			let ambient_light = 0.00;
			let final_intensity = ambient_light + diffuse_intensity * shadow_multiplier;
			surface_color = hit_sphere.albedo_and_atmos_flag.xyz * final_intensity * light_color;
		}
		
		final_pixel_color = surface_color * total_transmittance + total_in_scattering;

		// City lights from GOL texture (emissive, not attenuated)
		if (hit_sphere.emissive_and_terrain_flag.w > 0.5) {
			let p_on_sphere_normalized = normalize(rec.p - hit_sphere.pos_and_radius.xyz);
			let uv = clamp(worldPosToUV(p_on_sphere_normalized), vec2<f32>(0.0), vec2<f32>(1.0));
			let dims = vec2<f32>(textureDimensions(gol_texture));
			let coords = vec2<u32>(uv * dims);
			let cell_state = textureLoad(gol_texture, coords, 0).r;
			if (cell_state == 1u) {
				let city_color = vec3<f32>(1.0, 0.85, 0.6);
				let city_intensity = 1.5;
				final_pixel_color += city_color * city_intensity;
			}
		}
		final_alpha = 1.0;
	} else {
		final_pixel_color = total_in_scattering;
	}

	let hit_dist = select(1.0e10, rec.t, rec.hit);
	textureStore(depthOut, coords_i, vec4<f32>(hit_dist, 0.0, 0.0, 0.0));
	textureStore(output, coords_i, vec4<f32>(final_pixel_color, final_alpha));
}

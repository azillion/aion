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

#include "raytracing.wgsl"
#include "atmosphere.wgsl"
#include "shading.wgsl"
 

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

struct ComputeParams { bodyCount: u32 };
@group(0) @binding(0) var<uniform> params: ComputeParams;
@group(0) @binding(1) var<storage, read> spheres: array<Sphere>;
@group(0) @binding(2) var output: texture_storage_2d<rgba16float, write>;
@group(0) @binding(3) var<uniform> camera: CameraUniforms;
@group(0) @binding(4) var depthOut: texture_storage_2d<r32float, write>;
@group(0) @binding(5) var<uniform> scene: SceneUniforms;


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

    // Ray-scene intersection against the tier's spheres buffer
    var rec = hit_tier_spheres(ray, createInterval(0.001, INFINITY), &spheres, params.bodyCount, vec3<f32>(1e30));

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
                let self_pos_world_relative = atmos_sphere.pos_and_radius.xyz * scene.tier_scale_and_pad.x;
				let atmos = get_sky_color(
                    Ray(ray.origin - atmos_sphere.pos_and_radius.xyz, ray.direction),
                    atmos_sphere.pos_and_radius.w, atmos_sphere.terrain_params.x,
                    atmos_sphere.pos_and_radius.xyz,
					-normalize(scene.dominant_light_direction.xyz),
					scene.dominant_light_color_and_debug.xyz,
                    scene, max_dist, &shadow_casters, shadowParams, self_pos_world_relative
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
			let ground_color = shade_planet_surface(rec, hit_sphere, ray, scene, camera);
			let self_pos_world_relative = hit_sphere.pos_and_radius.xyz * scene.tier_scale_and_pad.x;
			let view_atmos = get_sky_color(
				Ray(ray.origin - hit_sphere.pos_and_radius.xyz, ray.direction),
				hit_sphere.pos_and_radius.w, hit_sphere.terrain_params.x,
				hit_sphere.pos_and_radius.xyz,
				-normalize(scene.dominant_light_direction.xyz),
				scene.dominant_light_color_and_debug.xyz,
				scene, rec.t, &shadow_casters, shadowParams, self_pos_world_relative
			);
			surface_color = ground_color * view_atmos.transmittance + view_atmos.in_scattering;
		} else { // Moon or other simple sphere
			let light_dir_to_source = -normalize(scene.dominant_light_direction.xyz);
			let light_color = scene.dominant_light_color_and_debug.xyz;
			let diffuse_intensity = max(dot(rec.normal, light_dir_to_source), 0.0);
			let p_world = rec.p * scene.tier_scale_and_pad.x;
            let shadow_ray = Ray(p_world + rec.normal * 0.2, light_dir_to_source);
            let self_pos_world_unscaled = hit_sphere.pos_and_radius.xyz * scene.tier_scale_and_pad.x;
            let shadow_rec = hit_spheres_shadow(shadow_ray, createInterval(0.001, INFINITY), &shadow_casters, shadowParams.count, self_pos_world_unscaled);
			var shadow_multiplier = 1.0;
            if (shadow_rec.hit && dot(shadow_rec.emissive, shadow_rec.emissive) < 0.1) { shadow_multiplier = 0.01; }
			let ambient_light = 0.00;
			let final_intensity = ambient_light + diffuse_intensity * shadow_multiplier;
			surface_color = hit_sphere.albedo_and_atmos_flag.xyz * final_intensity * light_color;
		}
		
		final_pixel_color = surface_color * total_transmittance + total_in_scattering;
		final_alpha = 1.0;
	} else {
		final_pixel_color = total_in_scattering;
	}

	let hit_dist = select(1.0e10, rec.t, rec.hit);
	textureStore(depthOut, coords_i, vec4<f32>(hit_dist, 0.0, 0.0, 0.0));
	textureStore(output, coords_i, vec4<f32>(final_pixel_color, final_alpha));
}
#include "camera.wgsl"
#include "sceneUniforms.wgsl"
#include "noise.wgsl"
#include "planetSdf.wgsl"
#include "atmosphere.wgsl"

const INFINITY: f32 = 1e38;

fn degreesToRadians(degrees: f32) -> f32 { return degrees * PI / 180.0; }

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

fn rayAt(ray: Ray, t: f32) -> vec3<f32> { return ray.origin + ray.direction * t; }

fn hit_sphere(sphere: Sphere, r: Ray, ray_t: Interval, index: u32) -> HitRecord {
    var rec: HitRecord; rec.hit = false;
    let center = sphere.pos_and_radius.xyz; let radius = sphere.pos_and_radius.w;
    let oc = center - r.origin; let a = dot(r.direction, r.direction); let h = dot(r.direction, oc); let c = dot(oc, oc) - radius * radius;
    let discriminant = h * h - a * c; if (discriminant < 0.0) { return rec; }
    let sqrtd = sqrt(discriminant);
    var root = (h - sqrtd) / a; if (!intervalSurrounds(ray_t, root)) { root = (h + sqrtd) / a; if (!intervalSurrounds(ray_t, root)) { return rec; } }
    rec.t = root; rec.p = rayAt(r, rec.t);
    let outward_normal = (rec.p - center) / radius;
    rec.front_face = dot(r.direction, outward_normal) < 0.0;
    rec.normal = outward_normal;
    rec.albedo = sphere.albedo_and_atmos_flag.xyz; rec.emissive = sphere.emissive_and_terrain_flag.xyz; rec.object_index = index; rec.hit = true; return rec;
}

fn hit_spheres(r: Ray, ray_t: Interval, ignore_index: u32) -> HitRecord {
    var closest_so_far = ray_t.maxI; var rec: HitRecord; rec.hit = false;
    for (var i = 0u; i < params.bodyCount; i++) {
        if (i == ignore_index) { continue; }
        let sphere_rec = hit_sphere(spheres[i], r, Interval(ray_t.minI, closest_so_far), i);
        if (sphere_rec.hit) { closest_so_far = sphere_rec.t; rec = sphere_rec; }
    }
    return rec;
}

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
    let terrain_uniforms = TerrainUniforms( sphere.terrain_params.x, sphere.terrain_params.y, sphere.terrain_params.z, sphere.terrain_params.w );
    let tier_scale = scene.tier_scale_and_pad.x;
    let R_world = terrain_uniforms.base_radius * tier_scale;
    let p_local_hit = rec.p - sphere.pos_and_radius.xyz;
    let h_world = h_noise(normalize(p_local_hit), terrain_uniforms, rec.t * tier_scale, R_world, camera);
    let is_ocean = h_world < terrain_uniforms.sea_level;

    var surface_albedo: vec3<f32>; var surface_normal: vec3<f32>;
    if (is_ocean) {
        surface_albedo = get_ocean_material(p_local_hit).albedo;
        surface_normal = get_ocean_wave_normal(p_local_hit, normalize(p_local_hit));
    } else {
        surface_normal = get_terrain_normal_from_heightfield( rec.p, sphere.pos_and_radius.xyz, rec.normal, terrain_uniforms, rec.t, R_world, camera, scene );
        if (terrain_uniforms.seed == 1337.0) {
            surface_albedo = get_lunar_material(p_local_hit, surface_normal).albedo;
        } else {
            surface_albedo = get_material(p_local_hit, surface_normal, terrain_uniforms, h_world).albedo;
        }
    }

    let light_dir_to_source = -normalize(scene.dominant_light_direction.xyz);
    let view_dir = normalize(ray.origin - rec.p);
    var shadow_multiplier = 1.0;
    let shadow_ray = Ray(rec.p + surface_normal * 0.02, light_dir_to_source);
    let inter_body_shadow_rec = hit_spheres(shadow_ray, createInterval(0.001, INFINITY), rec.object_index);
    if (inter_body_shadow_rec.hit && dot(inter_body_shadow_rec.emissive, inter_body_shadow_rec.emissive) < 0.1) { shadow_multiplier = 0.0; }

    let diffuse = max(0.0, dot(surface_normal, light_dir_to_source));
    let half_vec = normalize(light_dir_to_source + view_dir);
    let specular = pow(max(0.0, dot(surface_normal, half_vec)), 64.0) * select(0.1, 1.0, is_ocean);
    // Ambient tied to atmosphere presence (albedo_and_atmos_flag.w)
    let has_atmosphere = sphere.albedo_and_atmos_flag.w > 0.5;
    let ambient = select(0.0, 0.00, has_atmosphere);
    let lit_color = (surface_albedo * (diffuse + ambient) + specular) * shadow_multiplier;
    return lit_color * scene.dominant_light_color_and_debug.xyz;
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
	var rec = hit_spheres(ray, createInterval(0.001, INFINITY), 9999u);
	var final_pixel_color = vec3<f32>(0.0);
	var final_alpha = 0.0;

	if (rec.hit) {
		let hit_sphere = spheres[rec.object_index];
		var surface_color: vec3<f32>;
		if (dot(hit_sphere.emissive_and_terrain_flag.xyz, hit_sphere.emissive_and_terrain_flag.xyz) > 0.1) {
			surface_color = hit_sphere.emissive_and_terrain_flag.xyz;
		} else if (hit_sphere.emissive_and_terrain_flag.w > 0.5) {
			surface_color = shade_planet_surface(rec, hit_sphere, ray, scene, camera);
		} else {
			let light_dir_to_source = -normalize(scene.dominant_light_direction.xyz);
			let diffuse_intensity = max(dot(rec.normal, light_dir_to_source), 0.0);
			let shadow_ray = Ray(rec.p + rec.normal * 0.02, light_dir_to_source);
			let shadow_rec = hit_spheres(shadow_ray, createInterval(0.001, INFINITY), rec.object_index);
			var shadow_multiplier = 1.0;
			if (shadow_rec.hit && dot(shadow_rec.emissive, shadow_rec.emissive) < 0.1) { shadow_multiplier = 0.0; }
			let ambient_light = 0.0;
			let final_intensity = ambient_light + diffuse_intensity * shadow_multiplier;
			surface_color = hit_sphere.albedo_and_atmos_flag.xyz * final_intensity * scene.dominant_light_color_and_debug.xyz;
		}
		var total_in_scattering = vec3<f32>(0.0);
		var total_transmittance = vec3<f32>(1.0);
		if (atmosphereEnabled) {
			for (var i = 0u; i < params.bodyCount; i = i + 1u) {
				let atmos_sphere = spheres[i];
				if (atmos_sphere.albedo_and_atmos_flag.w > 0.5) {
					let atmos = get_sky_color(ray, atmos_sphere.pos_and_radius.xyz, atmos_sphere.terrain_params.x, -normalize(scene.dominant_light_direction.xyz), scene, rec.t);
					total_in_scattering += atmos.in_scattering;
					total_transmittance *= atmos.transmittance;
				}
			}
		}
		final_pixel_color = surface_color * total_transmittance + total_in_scattering;
		final_alpha = 1.0;
    } else {
        if (atmosphereEnabled) {
            for (var i = 0u; i < params.bodyCount; i = i + 1u) {
                let atmos_sphere = spheres[i];
			if (atmos_sphere.albedo_and_atmos_flag.w > 0.5) {
                    // Occlusion: skip if the solid planet blocks the view to its far-side atmosphere
                    let t_solid = ray_sphere_intersect(ray.origin - atmos_sphere.pos_and_radius.xyz, ray.direction, atmos_sphere.pos_and_radius.w);
                    if (t_solid.x > 0.0) { continue; }

                    let atmos = get_sky_color(ray, atmos_sphere.pos_and_radius.xyz, atmos_sphere.terrain_params.x, -normalize(scene.dominant_light_direction.xyz), scene, INFINITY);
                    final_pixel_color += atmos.in_scattering;
                    final_alpha = max(final_alpha, atmos.alpha);
                }
            }
        }
    }
	let hit_dist = select(1.0e10, rec.t, rec.hit);
	textureStore(depthOut, coords_i, vec4<f32>(hit_dist, 0.0, 0.0, 0.0));
	textureStore(output, coords_i, vec4<f32>(final_pixel_color, final_alpha));
}
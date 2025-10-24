#include "camera.wgsl"
#include "sceneUniforms.wgsl"

const PI: f32 = 3.1415926535897932385;
const INFINITY: f32 = 1e38;
const SEED: vec2<f32> = vec2<f32>(69.68, 4.20);
const MAX_DEPTH: u32 = 100;

fn lerp(a: vec3<f32>, b: vec3<f32>, t: f32) -> vec3<f32> {
    return a * (1.0 - t) + b * t;
}

fn degreesToRadians(degrees: f32) -> f32 { return degrees * PI / 180.0; }

fn hash(seed: vec2<u32>) -> u32 { var state = seed.x; state = state ^ (state << 13u); state = state ^ (state >> 17u); state = state ^ (state << 5u); state = state * 1597334677u; state = state ^ seed.y; state = state * 1597334677u; return state; }
fn rand(seed: vec2<u32>) -> f32 { return f32(hash(seed)) / 4294967295.0; }
fn randMinMax(seed: vec2<u32>, min: f32, max: f32) -> f32 { return min + (max - min) * rand(seed); }
fn randVec3(seed: vec2<u32>) -> vec3<f32> { return vec3<f32>(rand(seed), rand(seed + vec2<u32>(1u, 0u)), rand(seed + vec2<u32>(0u, 1u))); }

struct Interval { minI: f32, maxI: f32 };
fn createInterval(min: f32, max: f32) -> Interval { return Interval(min, max); }
fn intervalSurrounds(i: Interval, x: f32) -> bool { return i.minI < x && x < i.maxI; }
const INTERVAL_UNIVERSE: Interval = Interval(-INFINITY, INFINITY);

struct Sphere {
    pos_and_radius: vec4<f32>,
    albedo_and_pad: vec4<f32>,
    emissive_and_fuzz: vec4<f32>,
    ref_idx_and_opacity: vec4<f32>,
}

struct HitRecord {
    p: vec3<f32>,
    normal: vec3<f32>,
    t: f32,
    hit: bool,
    front_face: bool,
    albedo: vec3<f32>,
    emissive: vec3<f32>,
    object_index: u32,
}

struct Ray { origin: vec3<f32>, direction: vec3<f32> }

fn hit_sphere(sphere: Sphere, r: Ray, ray_t: Interval, index: u32) -> HitRecord {
    var rec: HitRecord; rec.hit = false;
    let center = sphere.pos_and_radius.xyz; let radius = sphere.pos_and_radius.w;
    let oc = center - r.origin; let a = dot(r.direction, r.direction); let h = dot(r.direction, oc); let c = dot(oc, oc) - radius * radius;
    let discriminant = h * h - a * c; if (discriminant < 0.0) { return rec; }
    let sqrtd = sqrt(discriminant);
    var root = (h - sqrtd) / a; if (!intervalSurrounds(ray_t, root)) { root = (h + sqrtd) / a; if (!intervalSurrounds(ray_t, root)) { return rec; } }
    rec.t = root; rec.p = rayAt(r, rec.t);
    let outward_normal = (rec.p - center) / radius; rec.front_face = dot(r.direction, outward_normal) < 0.0; rec.normal = select(-outward_normal, outward_normal, rec.front_face);
    rec.albedo = sphere.albedo_and_pad.xyz; rec.emissive = sphere.emissive_and_fuzz.xyz; rec.object_index = index; rec.hit = true; return rec;
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

fn rayAt(ray: Ray, t: f32) -> vec3<f32> { return ray.origin + ray.direction * t; }

struct Camera {
    origin: vec3<f32>, lower_left_corner: vec3<f32>, horizontal: vec3<f32>, vertical: vec3<f32>, samples_per_pixel: u32, vfov: f32, u: vec3<f32>, v: vec3<f32>, w: vec3<f32>,
}

fn createCamera(aspect_ratio: f32) -> Camera {
    let samples_per_pixel: u32 = 1u; let vfov = 25.0; let lookfrom = vec3<f32>(0.0, 0.0, 0.0); let vup = camera.up; let focus_distance = max(0.001, camera.distance_to_target);
    let theta = degreesToRadians(vfov); let h = tan(theta / 2.0); let viewport_height = 2.0 * h; let viewport_width = aspect_ratio * viewport_height;
    let w = -camera.forward; let u = normalize(cross(vup, w)); let v = cross(w, u);
    let origin = lookfrom; let horizontal = viewport_width * u; let vertical = viewport_height * v; let lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - w;
    return Camera(origin, lower_left_corner, horizontal, vertical, samples_per_pixel, vfov, u, v, w);
}

fn getRay(camera: Camera, s: f32, t: f32, _seed: vec2<u32>) -> Ray {
    return Ray(camera.origin, normalize(camera.lower_left_corner + s*camera.horizontal + t*camera.vertical - camera.origin));
}

fn rayColor(initial_ray: Ray) -> vec3<f32> {
    let rec = hit_spheres(initial_ray, createInterval(0.001, INFINITY), 9999u);
    if (!rec.hit) { return vec3<f32>(0.0, 0.0, 0.0); }
    if (dot(rec.emissive, rec.emissive) > 0.1) { return rec.emissive; }
    // Use stable precomputed direction (negated to point from surface to light)
    let light_dir = -scene.dominant_light_direction.xyz;
    let diffuse_intensity = max(dot(rec.normal, light_dir), 0.0);
    return rec.albedo * diffuse_intensity * scene.dominant_light_color_and_debug.xyz;
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
    let coords = vec2<u32>(id.xy);
    if (coords.x >= dims.x || coords.y >= dims.y) { return; }
    let aspect_ratio = f32(dims.x) / f32(dims.y);
    let cam = createCamera(aspect_ratio);

    let u = (f32(coords.x) + 0.5) / f32(dims.x);
    let v = 1.0 - (f32(coords.y) + 0.5) / f32(dims.y);
    let seed = vec2<u32>(coords.x + dims.x * coords.y, 0u);
    let ray = getRay(cam, u, v, seed);

    // Trace once and get full hit
    let rec = hit_spheres(ray, createInterval(0.001, INFINITY), 9999u);

    // Derive color: restore lighting with simple directional shadow
    var pixel_color = vec3<f32>(0.0);
    if (rec.hit) {
        if (dot(rec.emissive, rec.emissive) > 0.1) {
            pixel_color = rec.emissive;
        } else {
            // Use stable precomputed direction (negated to point from surface to light)
            let light_dir = -scene.dominant_light_direction.xyz;
            let diffuse_intensity = max(dot(rec.normal, light_dir), 0.0);

            // With scaled tiers, a small epsilon avoids shadow acne without artifacts
            const SHADOW_BIAS: f32 = 0.01;
            let shadow_ray = Ray(rec.p + rec.normal * SHADOW_BIAS, light_dir);
            let dist_to_light = 1.0e12; // effectively infinity for directional light
            let shadow_rec = hit_spheres(shadow_ray, createInterval(0.001, dist_to_light), 9999u);

            var shadow_multiplier = 1.0;
            if (shadow_rec.hit && dot(shadow_rec.emissive, shadow_rec.emissive) < 0.1) {
                shadow_multiplier = 0.0;
            }

            pixel_color = rec.albedo * diffuse_intensity * scene.dominant_light_color_and_debug.xyz * shadow_multiplier;
        }
    }

    let hit_dist = select(0.0, rec.t, rec.hit);
    textureStore(depthOut, vec2<i32>(coords), vec4<f32>(hit_dist, 0.0, 0.0, 0.0));
    textureStore(output, vec2<i32>(coords), vec4<f32>(pixel_color, 1.0));
}



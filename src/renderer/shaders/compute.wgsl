const PI: f32 = 3.1415926535897932385;
const INFINITY: f32 = 1e38;
const SEED: vec2<f32> = vec2<f32>(69.68, 4.20);
const MAX_DEPTH: u32 = 100;

// CameraUniforms is provided by a shared include (camera.wgsl)

fn lerp(a: vec3<f32>, b: vec3<f32>, t: f32) -> vec3<f32> {
    return a * (1.0 - t) + b * t;
}

fn degreesToRadians(degrees: f32) -> f32 {
    return degrees * PI / 180.0;
}

fn hash(seed: vec2<u32>) -> u32 {
    var state = seed.x;
    state = state ^ (state << 13u);
    state = state ^ (state >> 17u);
    state = state ^ (state << 5u);
    state = state * 1597334677u;
    state = state ^ seed.y;
    state = state * 1597334677u;
    return state;
}

fn rand(seed: vec2<u32>) -> f32 {
    return f32(hash(seed)) / 4294967295.0;
}

fn randMinMax(seed: vec2<u32>, min: f32, max: f32) -> f32 {
    return min + (max - min) * rand(seed);
}

fn randVec3(seed: vec2<u32>) -> vec3<f32> {
    return vec3<f32>(rand(seed), rand(seed + vec2<u32>(1u, 0u)), rand(seed + vec2<u32>(0u, 1u)));
}

fn randVec3MinMax(seed: vec2<u32>, min: f32, max: f32) -> vec3<f32> {
    return vec3<f32>(randMinMax(seed, min, max), randMinMax(seed + vec2<u32>(1u, 0u), min, max), randMinMax(seed + vec2<u32>(0u, 1u), min, max));
}

fn randInUnitDisk(seed: vec2<u32>) -> vec3<f32> {
    var tempseed = seed;
    loop {
        let p = vec3<f32>(randMinMax(tempseed, -1.0, 1.0), randMinMax(tempseed + vec2<u32>(1u, 0u), -1.0, 1.0), 0.0);
        if dot(p, p) < 1.0 {
            return p;
        }
        tempseed = vec2<u32>(hash(tempseed), hash(tempseed + vec2<u32>(1u, 1u)));
    }
}

struct Interval {
    minI: f32,
    maxI: f32,
}

fn createInterval(min: f32, max: f32) -> Interval {
    return Interval(min, max);
}

fn intervalSize(i: Interval) -> f32 {
    return i.maxI - i.minI;
}

fn intervalContains(i: Interval, x: f32) -> bool {
    return i.minI <= x && x <= i.maxI;
}

fn intervalSurrounds(i: Interval, x: f32) -> bool {
    return i.minI < x && x < i.maxI;
}

fn clampInterval(i: Interval, minI: f32, maxI: f32) -> f32 {
    return min(max(i.minI, minI), maxI);
}

const INTERVAL_EMPTY: Interval = Interval(INFINITY, -INFINITY);
const INTERVAL_UNIVERSE: Interval = Interval(-INFINITY, INFINITY);

// Flattened Sphere layout to match CPU-packed 16-float buffer
struct Sphere {
    center: vec3<f32>,
    radius: f32,
    albedo: vec3<f32>,
    _pad1: f32,
    emissive: vec3<f32>,
    fuzziness: f32,
    refraction_index: f32,
    mat_type: u32,
    _pad2: u32,
    _pad3: u32,
}

struct HitRecord {
    p: vec3<f32>,
    normal: vec3<f32>,
    t: f32,
    hit: bool,
    front_face: bool,
    // Flattened material properties
    albedo: vec3<f32>,
    emissive: vec3<f32>,
}

fn hit_sphere(sphere: Sphere, r: Ray, ray_t: Interval) -> HitRecord {
    var rec: HitRecord;
    rec.hit = false;

    let oc = sphere.center - r.origin;
    let a = dot(r.direction, r.direction);
    let h = dot(r.direction, oc);
    let c = dot(oc, oc) - sphere.radius * sphere.radius;

    let discriminant = h * h - a * c;
    if (discriminant < 0.0) {
        return rec;
    }

    let sqrtd = sqrt(discriminant);

    var root = (h - sqrtd) / a;
    if (!intervalSurrounds(ray_t, root)) {
        root = (h + sqrtd) / a;
        if (!intervalSurrounds(ray_t, root)) {
            return rec;
        }
    }

	rec.t = root;
    rec.p = rayAt(r, rec.t);
    let outward_normal = (rec.p - sphere.center) / sphere.radius;
    rec.front_face = dot(r.direction, outward_normal) < 0.0;
    rec.normal = select(-outward_normal, outward_normal, rec.front_face);
    // Copy material properties directly
    rec.albedo = sphere.albedo;
    rec.emissive = sphere.emissive;
    rec.hit = true;

    return rec;
}

fn hit_spheres(r: Ray, world: array<Sphere, NUM_SPHERES>, ray_t: Interval) -> HitRecord {
    var closest_so_far = ray_t.maxI;
    var rec: HitRecord;
    rec.hit = false;

    for (var i = 0u; i < NUM_SPHERES; i++) { 
        let sphere_rec = hit_sphere(world[i], r, createInterval(ray_t.minI, closest_so_far));
        if (sphere_rec.hit) {
            closest_so_far = sphere_rec.t;
            rec = sphere_rec;
        }
    }

    return rec;
}

struct Camera {
    origin: vec3<f32>,
    lower_left_corner: vec3<f32>,
    horizontal: vec3<f32>,
    vertical: vec3<f32>,
    samples_per_pixel: u32,
    vfov: f32,
    lookfrom: vec3<f32>,
    lookat: vec3<f32>,
    vup: vec3<f32>,
    defocus_angle: f32,
    focus_distance: f32,
    u: vec3<f32>,
    v: vec3<f32>,
    w: vec3<f32>,
    defocus_disk_u: vec3<f32>,
    defocus_disk_v: vec3<f32>,
}

fn createCamera(aspect_ratio: f32) -> Camera {
    let samples_per_pixel: u32 = 1u; // reduce noise and debug visibility
    let vfov = 25.0;
    let lookfrom = camera.eye; // (0,0,0)
    let vup = camera.up;
    let defocus_angle = 0.0; // disable DOF to avoid over-blur
    let focus_distance = max(0.001, camera.distance_to_target);

    let theta = degreesToRadians(vfov);
    let h = tan(theta / 2.0);
    let viewport_height = 2.0 * h * focus_distance;
    let viewport_width = aspect_ratio * viewport_height;

    // Use provided forward vector for basis
    let w = -camera.forward;
    let u = normalize(cross(vup, w));
    let v = cross(w, u);

    let origin = lookfrom;
    // Scale by focus distance so the image plane sits at the correct depth
    let horizontal = focus_distance * viewport_width * u;
    let vertical = focus_distance * viewport_height * v;
    let lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - focus_distance * w;

    let defocus_radius = focus_distance * tan(degreesToRadians(defocus_angle / 2.0));
    let defocus_disk_u = u * defocus_radius;
    let defocus_disk_v = v * defocus_radius;

    return Camera(origin, lower_left_corner, horizontal, vertical, 
                  samples_per_pixel, vfov, lookfrom, (lookfrom + camera.forward), vup, 
                  defocus_angle, focus_distance, u, v, w, defocus_disk_u, defocus_disk_v);
}

fn getRay(camera: Camera, s: f32, t: f32, seed: vec2<u32>) -> Ray {
    var rd: vec3<f32> = vec3<f32>(0.0, 0.0, 0.0);
    if (camera.defocus_angle > 0.0) {
        let p = randInUnitDisk(seed);
        rd = (camera.defocus_disk_u * p.x + camera.defocus_disk_v * p.y);
    }
    let offset = camera.u * rd.x + camera.v * rd.y;
    return Ray(
        camera.origin + offset,
        camera.lower_left_corner + s*camera.horizontal + t*camera.vertical - camera.origin - offset
    );
}

struct Ray {
    origin: vec3<f32>,
    direction: vec3<f32>
}

const R = cos(PI / 4.0);

fn rayColor(initial_ray: Ray, world: array<Sphere, NUM_SPHERES>, seed: vec2<u32>) -> vec3<f32> {
    // Fixed near clip for primary rays to avoid front-surface clipping at close range
    let t_min = 0.0001;
    let rec = hit_spheres(initial_ray, world, createInterval(t_min, INFINITY));
    if (!rec.hit) {
        return vec3<f32>(0.0, 0.0, 0.0);
    }
    if (dot(rec.emissive, rec.emissive) > 0.0) {
        return rec.emissive;
    }
    let light_pos = spheres[0].center;
    let light_dir = normalize(light_pos - rec.p);
    let diffuse_intensity = max(dot(rec.normal, light_dir), 0.0);
    const SHADOW_BIAS: f32 = 0.0001;
    let shadow_ray = Ray(rec.p + rec.normal * SHADOW_BIAS, light_dir);
    let shadow_rec = hit_spheres(shadow_ray, world, createInterval(0.001, length(light_pos - rec.p)));
    // Path is clear if we hit nothing OR we hit an emissive (the light itself)
    if (!shadow_rec.hit || dot(shadow_rec.emissive, shadow_rec.emissive) > 0.0) {
        let light_brightness = 5.0;
        return rec.albedo * diffuse_intensity * light_brightness;
    } else {
        return rec.albedo * 0.05;
    }
}

fn rayAt(ray: Ray, t: f32) -> vec3<f32> {
    return ray.origin + ray.direction * t;
}

// Scene is now generated on CPU and uploaded via storage buffer

@group(0) @binding(0) var<storage, read> spheres: array<Sphere, NUM_SPHERES>;
@group(0) @binding(1) var output: texture_storage_2d<rgba16float, write>;
@group(0) @binding(2) var<uniform> camera: CameraUniforms;

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let dims = textureDimensions(output);
    let coords = vec2<u32>(id.xy);
    
    if (coords.x >= dims.x || coords.y >= dims.y) {
        return;
    }

    let aspect_ratio = f32(dims.x) / f32(dims.y);
    let cam = createCamera(aspect_ratio);
    // Scene data is pre-populated via storage buffer

    var pixel_color = vec3<f32>(0.0, 0.0, 0.0);
    let s_limit = cam.samples_per_pixel;
    for (var s = 0u; s < s_limit; s++) {
        let seed = vec2<u32>(coords.x + dims.x * coords.y, s);
        let u = (f32(coords.x) + rand(seed)) / f32(dims.x);
        let v = 1.0 - (f32(coords.y) + rand(seed + vec2<u32>(1u, 1u))) / f32(dims.y);
        let ray = getRay(cam, u, v, seed);
        pixel_color += rayColor(ray, spheres, seed);
    }
    pixel_color = pixel_color / f32(s_limit);

    textureStore(output, vec2<i32>(coords), vec4<f32>(pixel_color, 1.0));
}



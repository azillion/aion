const PI: f32 = 3.1415926535897932385;
const INFINITY: f32 = 1e38;
const SEED: vec2<f32> = vec2<f32>(69.68, 4.20);
const MAX_DEPTH: u32 = 100;
// NUM_SPHERES is injected dynamically by the host code

// Camera uniforms provided by the host
struct CameraUniforms {
    eye: vec3<f32>,
    look_at: vec3<f32>,
    up: vec3<f32>,
}

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

fn randInUnitSphere(seed: vec2<u32>) -> vec3<f32> {
    var tempseed = seed;
    loop {
        let p = randVec3MinMax(tempseed, -1.0, 1.0);
        if length(p) < 1.0 {
            return p;
        }
        tempseed = vec2<u32>(hash(tempseed), hash(tempseed + vec2<u32>(1u, 1u)));
    }
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

fn randUnitVector(seed: vec2<u32>) -> vec3<f32> {
    return normalize(randInUnitSphere(seed));
}

fn randomOnHemisphere(normal: vec3<f32>, seed: vec2<u32>) -> vec3<f32> {
    let on_unit_sphere = randUnitVector(seed);
    if (dot(on_unit_sphere, normal) > 0.0) {
        return on_unit_sphere;
    } else {
        return -on_unit_sphere;
    }
}

fn reflect(v: vec3<f32>, n: vec3<f32>) -> vec3<f32> {
    return v - 2.0 * dot(v, n) * n;
}

fn reflectance(cosine: f32, ref_idx: f32) -> f32 {
    var r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow(1.0 - cosine, 5.0);
}

fn refract(uv: vec3<f32>, n: vec3<f32>, etai_over_etat: f32) -> vec3<f32> {
    let cos_theta = min(dot(-uv, n), 1.0);
    let r_out_perp = etai_over_etat * (uv + cos_theta * n);
    let r_out_parallel = -sqrt(abs(1.0 - length(r_out_perp) * length(r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}

fn nearZero(v: vec3<f32>) -> bool {
    let s = 1e-8;
    return (v.x > -s && v.x < s) && (v.y > -s && v.y < s) && (v.z > -s && v.z < s);
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

struct Material {
    albedo: vec3<f32>,
    emissive: vec3<f32>,
    fuzziness: f32,
    refraction_index: f32,
    mat_type: u32,
}

struct ScatterRecord {
    scattered: Ray,
    attenuation: vec3<f32>,
    is_scattered: bool,
}

struct Sphere {
    center: vec3<f32>,
    radius: f32,
    material: Material,
}

struct HitRecord {
    p: vec3<f32>,
    normal: vec3<f32>,
    t: f32,
    hit: bool,
    front_face: bool,
    material: Material,
}

struct FaceNormalRecord {
    front_face: bool,
    normal: vec3<f32>,
}

fn setFaceNormal(rec: HitRecord, r: Ray, outwardNormal: vec3<f32>) -> FaceNormalRecord {
    var newRec: FaceNormalRecord;
    let frontFaceDirections = dot(r.direction, outwardNormal);
    if (frontFaceDirections < 0.0) {
        newRec.front_face = true;
        newRec.normal = outwardNormal;
    } else {
        newRec.front_face = false;
        newRec.normal = -outwardNormal;
    }
    return newRec;
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
    if (!intervalContains(ray_t, root)) {
        root = (h + sqrtd) / a;
        if (!intervalContains(ray_t, root)) {
            return rec;
        }
    }

    rec.t = root;
    rec.p = rayAt(r, rec.t);
    rec.normal = (rec.p - sphere.center) / sphere.radius;
    let faceNormalRec = setFaceNormal(rec, r, rec.normal);
    rec.front_face = faceNormalRec.front_face;
    rec.normal = faceNormalRec.normal;
    rec.material = sphere.material;
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
    let lookfrom = camera.eye;
    let lookat = camera.look_at;
    let vup = camera.up;
    let defocus_angle = 0.0; // disable DOF to avoid over-blur
    let focus_distance = length(lookfrom - lookat);

    let theta = degreesToRadians(vfov);
    let h = tan(theta / 2.0);
    let viewport_height = 2.0 * h * focus_distance;
    let viewport_width = aspect_ratio * viewport_height;

    let w = normalize(lookfrom - lookat);
    let u = normalize(cross(vup, w));
    let v = cross(w, u);

    let origin = lookfrom;
    let horizontal = viewport_width * u;
    let vertical = viewport_height * v;
    let lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - focus_distance * w;

    let defocus_radius = focus_distance * tan(degreesToRadians(defocus_angle / 2.0));
    let defocus_disk_u = u * defocus_radius;
    let defocus_disk_v = v * defocus_radius;

    return Camera(origin, lower_left_corner, horizontal, vertical, 
                  samples_per_pixel, vfov, lookfrom, lookat, vup, 
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

// Removed unused material constants now that scene materials are provided via storage buffer

const R = cos(PI / 4.0);

fn scatterLambertian(r: Ray, rec: HitRecord, material: Material, seed: vec2<u32>) -> ScatterRecord {
    var scatter_direction = rec.normal + randUnitVector(seed);
    if (nearZero(scatter_direction)) {
        scatter_direction = rec.normal;
    }
    let scattered = Ray(rec.p, scatter_direction);
    let attenuation = material.albedo;
    return ScatterRecord(scattered, attenuation, true);
}

fn scatterMetal(r: Ray, rec: HitRecord, material: Material, seed: vec2<u32>) -> ScatterRecord {
    let reflected = reflect(normalize(r.direction), rec.normal);
    let scattered = Ray(rec.p, reflected + material.fuzziness * randInUnitSphere(seed));
    let attenuation = material.albedo;
    let is_scattered = dot(scattered.direction, rec.normal) > 0.0;
    return ScatterRecord(scattered, attenuation, is_scattered);
}

fn scatterDielectric(r: Ray, rec: HitRecord, material: Material, seed: vec2<u32>) -> ScatterRecord {
    let attenuation = vec3<f32>(1.0, 1.0, 1.0);
    var refraction_ratio: f32;
    if (rec.front_face) {
        refraction_ratio = 1.0 / material.refraction_index;
    } else {
        refraction_ratio = material.refraction_index;
    }

    let unit_direction = normalize(r.direction);
    let cos_theta = min(dot(-unit_direction, rec.normal), 1.0);
    let sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    let cannot_refract = refraction_ratio * sin_theta > 1.0;
    var direction: vec3<f32>;
    if (cannot_refract || reflectance(cos_theta, refraction_ratio) > rand(seed)) {
        direction = reflect(unit_direction, rec.normal);
    } else {
        direction = refract(unit_direction, rec.normal, refraction_ratio);
    }

    return ScatterRecord(Ray(rec.p, direction), attenuation, true);
}

fn rayColor(initial_ray: Ray, world: array<Sphere, NUM_SPHERES>, seed: vec2<u32>) -> vec3<f32> {
    var accumulated_color = vec3<f32>(0.0, 0.0, 0.0);
    var attenuation = vec3<f32>(1.0, 1.0, 1.0);
    var ray = initial_ray;
    var current_seed = seed;

    for (var depth = 0u; depth < MAX_DEPTH; depth++) {
        let rec = hit_spheres(ray, world, createInterval(0.001, INFINITY));

        if (rec.hit) {
            current_seed = vec2<u32>(hash(current_seed), depth);

            // Add light from emissive materials
            accumulated_color += rec.material.emissive * attenuation;

            var scatterRec: ScatterRecord;
            if (rec.material.mat_type == 0u) { // Lambertian
                scatterRec = scatterLambertian(ray, rec, rec.material, current_seed);
            } else if (rec.material.mat_type == 1u) { // Metal
                scatterRec = scatterMetal(ray, rec, rec.material, current_seed);
            } else { // Dielectric
                scatterRec = scatterDielectric(ray, rec, rec.material, current_seed);
            }

            if (scatterRec.is_scattered) {
                attenuation *= scatterRec.attenuation;
                ray = scatterRec.scattered;
            } else {
                break; // Ray absorbed
            }

        } else {
            // Background is black in this scene
            break;
        }

        // Russian Roulette termination
        if (depth > 4u) {
            let p = max(attenuation.x, max(attenuation.y, attenuation.z));
            if (rand(current_seed) > p) {
                break;
            }
            attenuation = attenuation / max(p, 1e-3);
        }
    }
    
    return accumulated_color;
}

fn rayAt(ray: Ray, t: f32) -> vec3<f32> {
    return ray.origin + ray.direction * t;
}

// Scene is now generated on CPU and uploaded via storage buffer

@group(0) @binding(0) var<storage, read> spheres: array<Sphere, NUM_SPHERES>;
@group(0) @binding(1) var output: texture_storage_2d<rgba8unorm, write>;
@group(0) @binding(2) var<uniform> camera: CameraUniforms;

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let dims = textureDimensions(output);
    let coords = vec2<u32>(id.xy);
    
    if (coords.x >= dims.x || coords.y >= dims.y) {
        return;
    }

    let aspect_ratio = f32(dims.x) / f32(dims.y);
    let camera = createCamera(aspect_ratio);
    // Scene data is pre-populated via storage buffer

    var pixel_color = vec3<f32>(0.0, 0.0, 0.0);
    for (var s = 0u; s < camera.samples_per_pixel; s++) {
        let seed = vec2<u32>(coords.x + dims.x * coords.y, s);
        let u = (f32(coords.x) + rand(seed)) / f32(dims.x);
        let v = (f32(coords.y) + rand(seed + vec2<u32>(1u, 1u))) / f32(dims.y);
        let ray = getRay(camera, u, v, seed);
        pixel_color += rayColor(ray, spheres, seed);
    }
    
    pixel_color = sqrt(pixel_color / f32(camera.samples_per_pixel));

    textureStore(output, vec2<i32>(coords), vec4<f32>(pixel_color, 1.0));
}



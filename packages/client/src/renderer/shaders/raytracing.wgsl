fn hit_sphere(sphere: Sphere, oc: vec3<f32>, r: Ray, ray_t: Interval, index: u32) -> HitRecord {
    var rec: HitRecord; rec.hit = false;
    let radius = sphere.pos_high_and_radius.w;
    // --- Numerically Stable Ray-Sphere Intersection ---
    let t_ca = dot(oc, r.direction);
    let d_sq = dot(oc, oc) - t_ca * t_ca; // squared distance from center to ray
    let radius_sq = radius * radius;
    if (d_sq > radius_sq) { return rec; }
    let t_hc = sqrt(radius_sq - d_sq);
    var root = t_ca - t_hc;
    if (!intervalSurrounds(ray_t, root)) {
        root = t_ca + t_hc;
        if (!intervalSurrounds(ray_t, root)) { return rec; }
    }
    rec.t = root; rec.p = rayAt(r, rec.t);
    let outward_normal = (rec.p - oc) / radius; // rec.p is camera-relative, oc is camera-relative center
    rec.front_face = dot(r.direction, outward_normal) < 0.0;
    rec.normal = outward_normal;
    rec.albedo = sphere.albedo_and_atmos_flag.xyz; rec.emissive = sphere.emissive_and_terrain_flag.xyz; rec.object_index = index; rec.hit = true; return rec;
}

// Performs a numerically stable check against a list of spheres and returns the closest hit.
// The incoming ray `r` is in camera-relative space.
fn hit_spheres_shadow(r: Ray, ray_t: Interval, spheres_in: ptr<storage, array<Sphere>, read>, sphere_count: u32, ignore_pos: vec3<f32>, camera: CameraUniforms) -> HitRecord {
    var closest_so_far = ray_t.maxI;
    var rec: HitRecord;
    rec.hit = false;
    for (var i = 0u; i < sphere_count; i++) {
        let sphere = (*spheres_in)[i];
        let sphere_pos_relative = sub_f64_to_f32(
            sphere.pos_high_and_radius.xyz, sphere.pos_low_and_unused.xyz,
            camera.eye_pos_high, camera.eye_pos_low
        );
        if (distance(sphere_pos_relative, ignore_pos) < 0.1) { continue; }
        let oc = sphere_pos_relative - r.origin;
        let sphere_rec = hit_sphere(sphere, oc, r, Interval(ray_t.minI, closest_so_far), i);
        if (sphere_rec.hit) {
            closest_so_far = sphere_rec.t;
            rec = sphere_rec;
        }
    }
    return rec;
}

// Generic indexed hit against the current tier's spheres buffer. Ray origin MUST be (0,0,0).
fn hit_scene_spheres(r: Ray, ray_t: Interval, spheres_in: ptr<storage, array<Sphere>, read>, sphere_count: u32, camera: CameraUniforms) -> HitRecord {
    var closest_so_far = ray_t.maxI;
    var rec: HitRecord; rec.hit = false;
    for (var i = 0u; i < sphere_count; i++) {
        let s = (*spheres_in)[i];
        // Calculate camera-relative position using f64 emulation.
        // The ray origin r.origin is (0,0,0) in this space, so oc = sphere_pos_relative.
        let sphere_pos_relative = sub_f64_to_f32(
            s.pos_high_and_radius.xyz, s.pos_low_and_unused.xyz,
            camera.eye_pos_high, camera.eye_pos_low
        );
        let sphere_rec = hit_sphere(s, sphere_pos_relative, r, Interval(ray_t.minI, closest_so_far), i);
        if (sphere_rec.hit) { closest_so_far = sphere_rec.t; rec = sphere_rec; }
    }
    return rec;
}


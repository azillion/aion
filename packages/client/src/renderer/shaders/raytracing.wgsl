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
fn hit_spheres_shadow(r: Ray, ray_t: Interval, spheres_in: ptr<storage, array<Sphere>, read>, sphere_count: u32, ignore_pos: vec3<f32>) -> HitRecord {
    var closest_so_far = ray_t.maxI;
    var rec: HitRecord;
    rec.hit = false;
    for (var i = 0u; i < sphere_count; i++) {
        let sphere = (*spheres_in)[i];
        if (distance(sphere.pos_and_radius.xyz, ignore_pos) < 0.1) { continue; }
        let oc = sphere.pos_and_radius.xyz - r.origin;
        let sphere_rec = hit_sphere(sphere, oc, r, Interval(ray_t.minI, closest_so_far), i);
        if (sphere_rec.hit) {
            closest_so_far = sphere_rec.t;
            rec = sphere_rec;
        }
    }
    return rec;
}

// Generic indexed hit against the current tier's spheres buffer
fn hit_tier_spheres(r: Ray, ray_t: Interval, spheres_in: ptr<storage, array<Sphere>, read>, sphere_count: u32, ignore_pos: vec3<f32>) -> HitRecord {
    var closest_so_far = ray_t.maxI;
    var rec: HitRecord; rec.hit = false;
    for (var i = 0u; i < sphere_count; i++) {
        let s = (*spheres_in)[i];
        if (distance(s.pos_and_radius.xyz, ignore_pos) < 0.1) { continue; }
        let oc = s.pos_and_radius.xyz - r.origin;
        let sphere_rec = hit_sphere(s, oc, r, Interval(ray_t.minI, closest_so_far), i);
        if (sphere_rec.hit) { closest_so_far = sphere_rec.t; rec = sphere_rec; }
    }
    return rec;
}


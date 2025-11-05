// Core ray-tracing intersection functions for analytic spheres.
// Generic, robust ray-sphere intersection that works for ANY ray origin.
fn hit_sphere_generic(sphere: Sphere, sphere_center: vec3<f32>, r: Ray, ray_t: Interval, index: u32) -> HitRecord {
    var rec: HitRecord; rec.hit = false;
    let radius = sphere.pos_high_and_radius.w;
    // Vector from ray origin to sphere center (origin-to-center)
    let oc = sphere_center - r.origin;
    // --- Numerically Stable Ray-Sphere Intersection ---
    let t_ca = dot(oc, r.direction);
    if (t_ca < 0.0) { return rec; }
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
    let outward_normal = (rec.p - sphere_center) / radius; // outward from center to hit point
    rec.front_face = dot(r.direction, outward_normal) < 0.0;
    rec.normal = outward_normal;
    rec.albedo = sphere.albedo_and_atmos_flag.xyz; rec.emissive = sphere.emissive_and_terrain_flag.xyz; rec.object_index = index; rec.hit = true; return rec;
}

// Fast-path intersection assuming ray origin is (0,0,0). Used for camera rays.
fn hit_sphere_camera_ray(sphere: Sphere, sphere_center: vec3<f32>, r: Ray, ray_t: Interval, index: u32) -> HitRecord {
    return hit_sphere_generic(sphere, sphere_center, r, ray_t, index);
}

// Performs a numerically stable check against a list of spheres and returns the closest hit.
// The incoming ray `r` is in camera-relative space.
fn hit_spheres_shadow(r: Ray, ray_t: Interval, spheres_in: ptr<storage, array<Sphere>, read>, sphere_count: u32, ignore_pos: vec3<f32>, camera: CameraUniforms) -> HitRecord {
    var closest_so_far = ray_t.maxI;
    var rec: HitRecord;
    rec.hit = false;
    for (var i = 0u; i < sphere_count; i = i + 1u) {
        let sphere = (*spheres_in)[i];
        let sphere_pos_relative = sub_f64_to_f32(
            sphere.pos_high_and_radius.xyz, sphere.pos_low_and_unused.xyz,
            camera.eye_pos_high, camera.eye_pos_low
        );
        if (distance(sphere_pos_relative, ignore_pos) < 0.1) { continue; }
        // Shadow rays have non-zero origins, so use generic function.
        let sphere_rec = hit_sphere_generic(sphere, sphere_pos_relative, r, Interval(ray_t.minI, closest_so_far), i);
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
    for (var i = 0u; i < sphere_count; i = i + 1u) {
        let s = (*spheres_in)[i];
        // Calculate camera-relative position using f64 emulation.
        let sphere_pos_relative = sub_f64_to_f32(
            s.pos_high_and_radius.xyz, s.pos_low_and_unused.xyz,
            camera.eye_pos_high, camera.eye_pos_low
        );
        // Camera rays originate at (0,0,0), so use fast path.
        let sphere_rec = hit_sphere_camera_ray(s, sphere_pos_relative, r, Interval(ray_t.minI, closest_so_far), i);
        if (sphere_rec.hit) { closest_so_far = sphere_rec.t; rec = sphere_rec; }
    }
    return rec;
}


// World intersection that supports both analytic spheres and SDF terrain via ray marching.
fn hit_world(r: Ray, ray_t: Interval, spheres_in: ptr<storage, array<Sphere>, read>, sphere_count: u32, camera: CameraUniforms, scene: SceneUniforms) -> HitRecord {
    var closest_so_far = ray_t.maxI;
    var rec: HitRecord; rec.hit = false;

    for (var i = 0u; i < sphere_count; i = i + 1u) {
        let s = (*spheres_in)[i];
        let sphere_pos_relative = sub_f64_to_f32(
            s.pos_high_and_radius.xyz, s.pos_low_and_unused.xyz,
            camera.eye_pos_high, camera.eye_pos_low
        );

        let has_terrain = s.emissive_and_terrain_flag.w > 0.5;
        var current_rec: HitRecord; current_rec.hit = false;

        if (has_terrain) {
            let params = TerrainUniforms(
                s.terrain_params.x, s.terrain_params.y,
                s.terrain_params.z, s.terrain_params.w
            );
            let dist_to_planet_center = length(sphere_pos_relative);

            // Construct high-precision local ray using f64 emulation
            let origin_x = sub_f64(vec2<f32>(camera.eye_pos_high.x, camera.eye_pos_low.x), vec2<f32>(s.pos_high_and_radius.x, s.pos_low_and_unused.x));
            let origin_y = sub_f64(vec2<f32>(camera.eye_pos_high.y, camera.eye_pos_low.y), vec2<f32>(s.pos_high_and_radius.y, s.pos_low_and_unused.y));
            let origin_z = sub_f64(vec2<f32>(camera.eye_pos_high.z, camera.eye_pos_low.z), vec2<f32>(s.pos_high_and_radius.z, s.pos_low_and_unused.z));
            let ray_local_f64 = Ray_f64(
                vec3<f32>(origin_x.x, origin_y.x, origin_z.x),
                vec3<f32>(origin_x.y, origin_y.y, origin_z.y),
                r.direction
            );

            current_rec = ray_march(
                ray_local_f64,
                closest_so_far,
                params,
                camera,
                scene,
                dist_to_planet_center
            );

            if (current_rec.hit) {
                // Convert local-space hit point back to camera-relative space
                current_rec.p = current_rec.p + sphere_pos_relative;
                current_rec.object_index = i;
                // Use base albedo; detailed material handled in shading
                current_rec.albedo = s.albedo_and_atmos_flag.xyz;
                current_rec.emissive = s.emissive_and_terrain_flag.xyz;
            }
        } else {
            current_rec = hit_sphere_camera_ray(s, sphere_pos_relative, r, Interval(ray_t.minI, closest_so_far), i);
        }

        if (current_rec.hit) {
            closest_so_far = current_rec.t;
            rec = current_rec;
        }
    }

    return rec;
}



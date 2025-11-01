fn shade_planet_surface(rec: HitRecord, sphere: Sphere, ray: Ray, scene: SceneUniforms, camera: CameraUniforms) -> vec3<f32> {
    let surface_albedo = sphere.albedo_and_atmos_flag.xyz;
    let surface_normal = rec.normal;
    let is_ocean = sphere.terrain_params.y > 0.0; // Simple ocean check based on sea level

    let light_dir_to_source = -normalize(scene.dominant_light_direction.xyz);
    let light_color = scene.dominant_light_color_and_debug.xyz;
    let view_dir = normalize(ray.origin - rec.p);
    let tier_scale = scene.tier_scale_and_pad.x;
    let p_world = rec.p * tier_scale;
    let self_pos_world = sphere.pos_and_radius.xyz * tier_scale;
    let shadow_ray = Ray(p_world + surface_normal * 0.2, light_dir_to_source);
    let inter_body_shadow_rec = hit_spheres_shadow(shadow_ray, createInterval(0.001, INFINITY), &shadow_casters, shadowParams.count, self_pos_world);
    let is_eclipsed = inter_body_shadow_rec.hit && dot(inter_body_shadow_rec.emissive, inter_body_shadow_rec.emissive) < 0.1;

    let has_atmosphere = sphere.albedo_and_atmos_flag.w > 0.5;
    var lit_color = surface_albedo * select(0.0, 0.02, has_atmosphere);

    if (has_atmosphere) {
        let zenith_dir = surface_normal;
        let sky_ambient = get_sky_color(
            Ray(rec.p - sphere.pos_and_radius.xyz, zenith_dir),
            sphere.pos_and_radius.w, sphere.terrain_params.x,
            sphere.pos_and_radius.xyz, light_dir_to_source, light_color, scene, INFINITY,
            &shadow_casters, shadowParams, self_pos_world
        );
        lit_color = surface_albedo * sky_ambient.in_scattering;

        let dot_nl = dot(surface_normal, light_dir_to_source);
        if (dot_nl > 0.0 && !is_eclipsed) {
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

             let diffuse_term = surface_albedo * dot_nl * final_sun_color;
             let half_vec = normalize(light_dir_to_source + view_dir);
             let specular_term = pow(max(0.0, dot(surface_normal, half_vec)), 64.0) * select(0.1, 1.0, is_ocean) * final_sun_color;
             lit_color += (diffuse_term + specular_term);
        }
    }
    
    return lit_color;
}


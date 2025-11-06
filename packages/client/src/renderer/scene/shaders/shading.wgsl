// Shading and lighting calculations for planetary surfaces.

fn get_ocean_wave_normal(p_local: vec3<f32>, base_normal: vec3<f32>) -> vec3<f32> {
    let wave_freq = 50.0;
    let wave_amp = 0.01;
    let wave_noise = snoise(p_local * wave_freq);
    return normalize(base_normal + wave_amp * vec3<f32>(wave_noise));
}

fn shade_planet_surface(rec: HitRecord, sphere: Sphere, sphere_pos_relative: vec3<f32>, ray: Ray, scene: SceneUniforms, camera: CameraUniforms, water_state: texture_2d_array<f32>) -> vec3<f32> {
    let R = sphere.terrain_params.x;
    let p_local = rec.p - sphere_pos_relative;
    let dir = normalize(p_local);

    // Recompute procedural terrain height at this location
    let terrain_params = TerrainUniforms(sphere.terrain_params.x, sphere.terrain_params.y, sphere.terrain_params.z, sphere.terrain_params.w);
    let terrain_height = h_noise(dir, terrain_params, rec.t, R, camera, scene);

    // Sample dynamic water level
    let water_uv = directionToCubeUV(dir);
    let dims = textureDimensions(water_state, 0);
    let tex_coords = min(vec2<i32>(water_uv.xy * vec2<f32>(dims)), vec2<i32>(dims) - 1);
    let water_depth = textureLoad(water_state, tex_coords, i32(water_uv.z), 0).r;
    let water_surface_height = sphere.terrain_params.y + water_depth;

    var surface_albedo: vec3<f32>;
    var surface_normal: vec3<f32>;
    let is_ocean = terrain_height < water_surface_height;

    if (is_ocean) {
        surface_albedo = vec3<f32>(0.05, 0.15, 0.2);
        surface_normal = get_ocean_wave_normal(p_local, rec.normal);
    } else {
        surface_albedo = sphere.albedo_and_atmos_flag.xyz;
        surface_normal = rec.normal;
    }

    let light_dir_to_source = normalize(scene.dominant_light_direction.xyz);
    let light_color = scene.dominant_light_color_and_debug.xyz;
    let view_dir = normalize(ray.origin - rec.p);
    let self_pos_world = sphere_pos_relative;
    let shadow_ray = Ray(rec.p + surface_normal * 0.2, light_dir_to_source);
    let inter_body_shadow_rec = hit_spheres_shadow(shadow_ray, createInterval(0.001, INFINITY), &shadow_casters, shadowParams.count, self_pos_world, camera);
    let is_eclipsed = inter_body_shadow_rec.hit && dot(inter_body_shadow_rec.emissive, inter_body_shadow_rec.emissive) < 0.1;

    let has_atmosphere = sphere.albedo_and_atmos_flag.w > 0.5 && scene.scale_and_flags.y > 0.5;
    var lit_color = surface_albedo * 0.05; // modest ambient term

    if (has_atmosphere) {
        let dot_nl = dot(surface_normal, light_dir_to_source);
        if (dot_nl > 0.0 && !is_eclipsed) {
             let params = get_earth_atmosphere();
             let p_start_local = p_local;
             let optical_lengths = calculate_optical_length_to_space(
                 p_start_local, light_dir_to_source,
                 sphere.pos_high_and_radius.w, sphere.pos_high_and_radius.w * ATMOSPHERE_RADIUS_SCALE,
                 params.h_rayleigh, params.h_mie, scene.scale_and_flags.x
             );
             let optical_depth = (params.beta_rayleigh * optical_lengths.x + params.beta_mie * optical_lengths.y);
             let transmittance = exp(-optical_depth);
             let final_sun_color = light_color * transmittance;

             let diffuse_term = surface_albedo * dot_nl * final_sun_color;
             let half_vec = normalize(light_dir_to_source + view_dir);
             let specular_term = pow(max(0.0, dot(surface_normal, half_vec)), 64.0) * select(0.1, 1.0, is_ocean) * final_sun_color;
             lit_color += (diffuse_term + specular_term);
        }
    } else {
        let dot_nl = max(0.0, dot(surface_normal, light_dir_to_source));
        var shadow_mult = 1.0;
        if (is_eclipsed) { shadow_mult = 0.0; }
        let diffuse = surface_albedo * dot_nl * shadow_mult * light_color;
        lit_color += diffuse;
    }
    
    return lit_color;
}



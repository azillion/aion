// Simple atmosphere utilities

fn get_sky_color(ray_dir: vec3<f32>) -> vec3<f32> {
    // Basic blue sky that fades with elevation; space fades to black
    let t = pow(max(ray_dir.y, 0.0), 0.5);
    return vec3<f32>(0.2, 0.4, 0.8) * t;
}

// Placeholder single-scattering approximation for in-scattering along the view ray to the hit
fn get_in_scattering(
    ray: Ray,
    hit_dist: f32,
    sun_dir: vec3<f32>,
    planet_pos: vec3<f32>,
    R: f32,
    R_atmos: f32
) -> vec3<f32> {
    // For now, return a simple haze proportional to distance, lightly clamped
    // This verifies wiring and prevents overbright results
    let haze = clamp(hit_dist / 5000.0, 0.0, 0.2);
    return vec3<f32>(0.2, 0.4, 0.8) * haze;
}



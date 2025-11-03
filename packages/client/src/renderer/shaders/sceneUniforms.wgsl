struct SceneUniforms {
    dominant_light_direction: vec4<f32>,
    dominant_light_color_and_debug: vec4<f32>,
    scale_and_flags: vec4<f32>, // .x = scale (always 1.0), .y = showAtmosphere, .z, .w unused
};



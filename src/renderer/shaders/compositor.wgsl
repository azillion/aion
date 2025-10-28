//#include support for Scene uniforms
#include "sceneUniforms.wgsl"
@group(0) @binding(0) var nearColor: texture_2d<f32>;
@group(0) @binding(1) var midColor: texture_2d<f32>;
@group(0) @binding(2) var farColor: texture_2d<f32>;
@group(0) @binding(3) var nearDepth: texture_2d<f32>;
@group(0) @binding(4) var midDepth: texture_2d<f32>;
@group(0) @binding(5) var farDepth: texture_2d<f32>;
@group(0) @binding(6) var s: sampler;
@group(0) @binding(7) var<uniform> scene: SceneUniforms;

@vertex
fn vertexMain(@builtin(vertex_index) vertexIndex: u32) -> @builtin(position) vec4<f32> {
    let pos = array<vec2<f32>, 3>(vec2(-1.0, -3.0), vec2(3.0, 1.0), vec2(-1.0, 1.0));
    return vec4<f32>(pos[vertexIndex], 0.0, 1.0);
}

@fragment
fn fragmentMain(@builtin(position) fragCoord: vec4<f32>) -> @location(0) vec4<f32> {
    let uv = fragCoord.xy / vec2<f32>(textureDimensions(nearColor));
    let coords_i = vec2<i32>(fragCoord.xy);

    // 1. Sample depth from all tiers using textureLoad (correct for unfilterable textures).
    let near_d = textureLoad(nearDepth, coords_i, 0).r;
    let mid_d = textureLoad(midDepth, coords_i, 0).r;
    let far_d = textureLoad(farDepth, coords_i, 0).r;

    // 2. Find the minimum (closest) depth value.
    let min_depth = min(min(near_d, mid_d), far_d);
    
    // 3. Select the color from the tier that contains the closest object.
    var final_color = vec4<f32>(0.0);
    if (min_depth < 1.0e9) { // Check if we hit any solid object at all
        // Use textureSampleLevel for the filterable color textures.
        if (near_d <= min_depth + 0.001) { // Add a small epsilon for float comparison
            final_color = textureSampleLevel(nearColor, s, uv, 0.0);
        } else if (mid_d <= min_depth + 0.001) {
            final_color = textureSampleLevel(midColor, s, uv, 0.0);
        } else {
            final_color = textureSampleLevel(farColor, s, uv, 0.0);
        }
    } else {
        // 4. If no solid object was hit, composite the skies (alpha blend back-to-front).
        let far_c = textureSampleLevel(farColor, s, uv, 0.0);
        let mid_c = textureSampleLevel(midColor, s, uv, 0.0);
        let near_c = textureSampleLevel(nearColor, s, uv, 0.0);

        // Blend far sky -> mid sky -> near sky
        final_color = far_c;
        final_color = mix(final_color, mid_c, mid_c.a);
        final_color = mix(final_color, near_c, near_c.a);
    }
    
    return vec4<f32>(final_color.rgb, 1.0);
}


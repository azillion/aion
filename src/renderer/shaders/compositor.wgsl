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

    // Perform a back-to-front composite to correctly layer skies over distant solids.
    let near_c = textureSampleLevel(nearColor, s, uv, 0.0);
    let mid_c = textureSampleLevel(midColor, s, uv, 0.0);
    let far_c = textureSampleLevel(farColor, s, uv, 0.0);

    // Start with the far tier as the base layer
    var final_color = far_c;

    // Mid tier compositing: occlude if closer than far, otherwise blend its sky
    if (mid_d < far_d) {
        final_color = mid_c;
    } else {
        final_color = mix(final_color, mid_c, mid_c.a);
    }

    // Near tier compositing: occlude if closest, otherwise blend its sky
    let current_closest_d = min(mid_d, far_d);
    if (near_d < current_closest_d) {
        final_color = near_c;
    } else {
        final_color = mix(final_color, near_c, near_c.a);
    }

    return vec4<f32>(final_color.rgb, 1.0);
}


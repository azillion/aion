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
    var final_color = vec4<f32>(0.0, 0.0, 0.0, 1.0);

    switch (i32(scene.dominant_light_color_and_debug.w)) {
        case 0: {
            final_color = textureSampleLevel(nearColor, s, uv, 0.0);
            break;
        }
        case 1: {
            final_color = textureSampleLevel(midColor, s, uv, 0.0);
            break;
        }
        case 2: {
            final_color = textureSampleLevel(farColor, s, uv, 0.0);
            break;
        }
        default: {
            let near_c = textureSampleLevel(nearColor, s, uv, 0.0);
            if (near_c.a > 0.0) {
                final_color = near_c;
            } else {
                let mid_c = textureSampleLevel(midColor, s, uv, 0.0);
                if (mid_c.a > 0.0) {
                    final_color = mid_c;
                } else {
                    let far_c = textureSampleLevel(farColor, s, uv, 0.0);
                    if (far_c.a > 0.0) {
                        final_color = far_c;
                    } else {
                        final_color = vec4<f32>(0.0, 0.0, 0.0, 1.0);
                    }
                }
            }
        }
    }

    return final_color;
}



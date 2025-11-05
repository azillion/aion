// One-time terrain baking pass: write procedural terrain height into a texture.
#include "camera.wgsl"
#include "sceneUniforms.wgsl"
#include "cubeSphere.wgsl"
#include "coarseGrid.wgsl"
#include "noise.wgsl"
#include "terrainCommon.wgsl"
#include "planetHeight.wgsl"

@group(0) @binding(0) var output: texture_storage_2d_array<r32float, write>;
@group(0) @binding(1) var<uniform> camera: CameraUniforms;
@group(0) @binding(2) var<uniform> scene: SceneUniforms;
@group(0) @binding(5) var<uniform> planet_params: TerrainUniforms;

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let dims = textureDimensions(output);
    if (id.x >= dims.x || id.y >= dims.y) { return; }

    let uv = (vec2f(id.xy) + 0.5) / vec2f(dims.xy);
    let dir = cubeUVToDirection(uv, id.z);

    let height = h_noise(dir, planet_params, 1.0, planet_params.base_radius, camera, scene);
    textureStore(output, vec2i(id.xy), i32(id.z), vec4f(height, 0.0, 0.0, 0.0));
}



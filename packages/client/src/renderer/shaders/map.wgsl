#include "camera.wgsl"
// CameraUniforms is provided by a shared include (camera.wgsl)
@group(0) @binding(0) var<uniform> camera: CameraUniforms;

// This struct must match the layout of the spheresBuffer in the compute pass.
struct Sphere {
	pos_and_radius: vec4<f32>,              // .xyz = position, .w = geometric radius
	albedo_and_atmos_flag: vec4<f32>,       // .xyz = albedo, .w = has_atmosphere_flag
	emissive_and_terrain_flag: vec4<f32>,   // .xyz = emissive, .w = has_terrain_flag
	ref_idx_opacity_pad: vec4<f32>,         // .x = ref_idx, .y = opacity, .zw unused (placeholder for alignment)
	terrain_params: vec4<f32>,              // .x = base_radius, .y = sea_level, .z = max_height, .w = seed
	padding_vec: vec4<f32>,                 // Padding to maintain 24-float stride
}
@group(0) @binding(1) var<storage, read> bodies: array<Sphere>;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) color: vec3<f32>,
    @location(1) uv: vec2<f32>,
};

@vertex
fn vertexMain(
    @builtin(vertex_index) vertexIndex: u32,
    @builtin(instance_index) instanceIndex: u32
) -> VertexOutput {
    let body = bodies[instanceIndex];

    // Define a camera-facing quad (billboard)
    let corners = array<vec2<f32>, 4>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>(-1.0,  1.0),
        vec2<f32>( 1.0,  1.0)
    );
    let corner = corners[vertexIndex % 4u];

    // Project center to clip space, then offset in NDC for screen-aligned quad (constant on-screen size)
    let p_clip = camera.viewProjection * vec4<f32>(body.pos_and_radius.xyz, 1.0);
    let size_ndc = 0.006; // tweak visual size as needed
    let offset = corner * size_ndc * p_clip.w;

    var output: VertexOutput;
    output.position = vec4<f32>(p_clip.x + offset.x, p_clip.y + offset.y, p_clip.z, p_clip.w);
    // Use the emissive color for suns, albedo for planets.
    let albedo = body.albedo_and_atmos_flag.xyz;
    let emissive = body.emissive_and_terrain_flag.xyz;
    output.color = select(albedo, emissive, dot(emissive, emissive) > 0.1);
    output.uv = corner; // Pass quad coordinates to fragment shader.
    return output;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    // Create a soft, circular shape inside the quad.
    let dist = length(in.uv);
    let alpha = 1.0 - smoothstep(0.8, 1.0, dist);

    if (alpha < 0.01) {
        discard;
    }

    return vec4<f32>(in.color, alpha);
}



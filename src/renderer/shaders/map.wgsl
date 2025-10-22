struct Camera {
    viewProjection: mat4x4<f32>,
    right: vec4<f32>,
    up: vec4<f32>,
    _pad: vec4<f32>,
};
@group(0) @binding(0) var<uniform> camera: Camera;

// This struct must match the layout of the spheresBuffer in the compute pass.
struct Sphere {
    center: vec3<f32>,
    radius: f32,
    albedo: vec3<f32>,
    _pad1: f32,
    emissive: vec3<f32>,
    _pad2: f32,
    fuzziness: f32,
    refraction_index: f32,
    mat_type: u32,
    _pad3: u32,
};
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
    let p_clip = camera.viewProjection * vec4<f32>(body.center, 1.0);
    let size_ndc = 0.006; // tweak visual size as needed
    let offset = corner * size_ndc * p_clip.w;

    var output: VertexOutput;
    output.position = vec4<f32>(p_clip.x + offset.x, p_clip.y + offset.y, p_clip.z, p_clip.w);
    // Use the emissive color for suns, albedo for planets.
    output.color = select(body.albedo, body.emissive, dot(body.emissive, body.emissive) > 0.1);
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



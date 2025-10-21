struct Camera {
    viewProjection: mat4x4<f32>,
    // We add dummy vec4s to ensure the struct size is a multiple of 256 bytes, a common requirement.
    _dummy1: vec4<f32>,
    _dummy2: vec4<f32>,
    _dummy3: vec4<f32>,
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
    
    // We'll scale the quad by a fixed amount for now to make planets visible on the map.
    let billboard_radius = 0.75;

    // We need the camera's right and up vectors to orient the quad.
    // For an orthographic top-down view, these are simple world axes.
    let right = vec3<f32>(1.0, 0.0, 0.0);
    let up = vec3<f32>(0.0, 1.0, 0.0);

    let worldPos = body.center + (right * corner.x * billboard_radius) + (up * corner.y * billboard_radius);
    
    var output: VertexOutput;
    output.position = camera.viewProjection * vec4<f32>(worldPos, 1.0);
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



// CameraUniforms is provided by a shared include (camera.wgsl)
@group(0) @binding(0) var<uniform> camera: CameraUniforms;

// Per-instance data for each star
// Packed to 32-byte stride (8 floats):
// position.xyz, (implicit padding), color.xyz, size
struct Star {
    position: vec3<f32>,
    color: vec3<f32>,
    size: f32,
};
@group(0) @binding(1) var<storage, read> stars: array<Star>;

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
    let star = stars[instanceIndex];
    let corners = array<vec2<f32>, 4>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 1.0, -1.0),
        vec2<f32>(-1.0,  1.0),
        vec2<f32>( 1.0,  1.0)
    );
    let corner = corners[vertexIndex % 4u];

    let worldPos = star.position + (camera.right * corner.x * star.size) + (camera.up * corner.y * star.size);
    
    var output: VertexOutput;
    output.position = camera.viewProjection * vec4<f32>(worldPos, 1.0);
    output.color = star.color;
    output.uv = corner;
    return output;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    // Create a soft, circular shape for the star
    let dist = length(in.uv);
    let alpha = 1.0 - smoothstep(0.4, 0.5, dist);

    if (alpha < 0.01) {
        discard;
    }

    return vec4<f32>(in.color * alpha, alpha);
}



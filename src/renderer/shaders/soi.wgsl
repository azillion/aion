// CameraUniforms is provided by a shared include (camera.wgsl)
@group(0) @binding(0) var<uniform> camera: CameraUniforms;

struct SoiInstance {
    position: vec3<f32>,
    radius: f32,
};
@group(0) @binding(1) var<storage, read> instances: array<SoiInstance>;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
};

@vertex
fn vertexMain(
    @builtin(vertex_index) vertexIndex: u32,
    @builtin(instance_index) instanceIndex: u32
) -> VertexOutput {
    let instance = instances[instanceIndex];
    let corners = array<vec2<f32>, 4>(
        vec2<f32>(-1.0, -1.0), vec2<f32>(1.0, -1.0), vec2<f32>(-1.0, 1.0), vec2<f32>(1.0, 1.0)
    );
    let corner = corners[vertexIndex % 4u];

    // Billboard in world space using camera right/up, scaled by world radius
    let right = camera.right;
    let up = camera.up;
    let worldPos = instance.position + (right * corner.x * instance.radius) + (up * corner.y * instance.radius);

    var out: VertexOutput;
    out.position = camera.viewProjection * vec4<f32>(worldPos, 1.0);
    out.uv = corner;
    return out;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let dist = length(in.uv);
    if (dist > 1.0) { discard; }

    let edge_glow = smoothstep(1.0, 0.95, dist);
    let fill = smoothstep(0.98, 0.97, dist) * 0.1;
    let alpha = max(edge_glow, fill) * 0.5;
    return vec4<f32>(0.8, 0.8, 1.0, alpha);
}



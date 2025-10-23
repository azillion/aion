// CameraUniforms is provided by a shared include (camera.wgsl)
@group(0) @binding(0) var<uniform> camera: CameraUniforms;

struct GlyphInstance {
    position: vec3<f32>,
    scale: f32,
};
@group(0) @binding(1) var<storage, read> instances: array<GlyphInstance>;

@vertex
fn vertexMain(
    @builtin(vertex_index) vertexIndex: u32,
    @builtin(instance_index) instanceIndex: u32
) -> @builtin(position) vec4<f32> {
    let instance = instances[instanceIndex];

    // Base triangle in NDC space (unit size), we'll scale to a small screen size
    let pos = array<vec2<f32>, 3>(
        vec2<f32>( 0.0,  0.5),
        vec2<f32>(-0.5, -0.5),
        vec2<f32>( 0.5, -0.5)
    );

    let p_clip = camera.viewProjection * vec4<f32>(instance.position, 1.0);
    let size = 0.004 * abs(instance.scale);
    let flip = select(-1.0, 1.0, instance.scale >= 0.0);
    let base = vec2<f32>(pos[vertexIndex].x, pos[vertexIndex].y * flip);
    let corner = base * size * p_clip.w;

    // Offset in clip space (similar to orbit ribbons)
    return vec4<f32>(p_clip.x + corner.x, p_clip.y + corner.y, p_clip.z, p_clip.w);
}

@fragment
fn fragmentMain() -> @location(0) vec4<f32> {
    return vec4<f32>(1.0, 1.0, 1.0, 1.0);
}



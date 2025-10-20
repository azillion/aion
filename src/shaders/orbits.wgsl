struct OrbitsUniforms {
    viewProjection: mat4x4<f32>,
    color: vec3<f32>,
    // implicit padding to 16 bytes
};

@group(0) @binding(0) var<uniform> uniforms: OrbitsUniforms;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
};

@vertex
fn vertexMain(@location(0) position: vec3<f32>) -> VertexOutput {
    var out: VertexOutput;
    out.position = uniforms.viewProjection * vec4<f32>(position, 1.0);
    return out;
}

@fragment
fn fragmentMain() -> @location(0) vec4<f32> {
    return vec4<f32>(uniforms.color, 1.0);
}



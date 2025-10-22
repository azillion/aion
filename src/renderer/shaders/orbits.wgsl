struct OrbitsUniforms {
    viewProjection: mat4x4<f32>,
    color: vec3<f32>,
    pointCount: f32,
    semiMajorAxis: f32,
    eccentricity: f32,
    currentTrueAnomaly: f32,
    // padding
    _pad1: f32,
    _pad2: f32,
};

@group(0) @binding(0) var<uniform> uniforms: OrbitsUniforms;
struct OrbitPointsBuffer {
    // 16-byte alignment for storage buffers
    points: array<vec4<f32>>,
};
@group(0) @binding(1) var<storage, read> orbit_points: OrbitPointsBuffer;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) point_index: f32,
};

@vertex
fn vertexMain(@builtin(vertex_index) vertexIndex: u32) -> VertexOutput {
    let point_index = vertexIndex / 2u;
    let is_right_side = (vertexIndex % 2u) == 1u;
    let count = u32(uniforms.pointCount);

    let i0 = select(point_index, point_index - 1u, point_index > 0u);
    let i2 = select(point_index + 1u, count - 1u, point_index + 1u >= count);

    let p_prev = orbit_points.points[i0].xyz;
    let p_curr = orbit_points.points[point_index].xyz;
    let p_next = orbit_points.points[i2].xyz;

    // Project to clip space
    let prev_clip = uniforms.viewProjection * vec4<f32>(p_prev, 1.0);
    let curr_clip = uniforms.viewProjection * vec4<f32>(p_curr, 1.0);
    let next_clip = uniforms.viewProjection * vec4<f32>(p_next, 1.0);

    if (curr_clip.w <= 0.0) {
        return VertexOutput(vec4<f32>(2.0, 2.0, 2.0, 1.0), 0.0);
    }

    let prev_ndc = prev_clip.xy / prev_clip.w;
    let curr_ndc = curr_clip.xy / curr_clip.w;
    let next_ndc = next_clip.xy / next_clip.w;

    var dir = normalize(next_ndc - prev_ndc);
    if (length(dir) == 0.0) {
        dir = vec2<f32>(1.0, 0.0);
    }

    let normal = vec2<f32>(-dir.y, dir.x);
    // Calculate variable ribbon width
    let a = uniforms.semiMajorAxis;
    let e = uniforms.eccentricity;
    let r_p = a * (1.0 - e); // Periapsis distance
    let r_a = a * (1.0 + e); // Apoapsis distance
    let r_curr = length(p_curr); // Current distance from barycenter

    // Normalize current distance to a 0..1 range between periapsis and apoapsis
    let t = (r_curr - r_p) / max(0.001, r_a - r_p);

    let maxWidth = 0.004;
    let minWidth = 0.001;
    // mix is thick->thin, so we want t=0 (periapsis) to be maxWidth
    let ribbon_width_ndc = mix(maxWidth, minWidth, t);
    let side_offset = select(-1.0, 1.0, is_right_side);

    var final_pos_clip = curr_clip;
    let offset = normal * ribbon_width_ndc * side_offset * final_pos_clip.w;
    final_pos_clip = vec4<f32>(final_pos_clip.x + offset.x, final_pos_clip.y + offset.y, final_pos_clip.z, final_pos_clip.w);

    return VertexOutput(final_pos_clip, f32(point_index));
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let total_points = uniforms.pointCount;
    let current_point_angle = (in.point_index / total_points) * 2.0 * 3.14159;

    // Normalize angles to be in 0..2PI range for comparison
    var current_nu_normalized = uniforms.currentTrueAnomaly % (2.0 * 3.14159);
    if (current_nu_normalized < 0.0) {
        current_nu_normalized += 2.0 * 3.14159;
    }

    let is_past = current_point_angle < current_nu_normalized;

    if (is_past) {
        // Discard fragments to create a dotted line effect for the past
        if (i32(in.point_index) % 8 < 4) {
            discard;
        }
    }

    return vec4<f32>(uniforms.color, 1.0);
}

// Post-processing and final compositing shader.
struct Theme {
    // Use vec4 for all fields to guarantee 16-byte alignment per field.
    // params.x = deltaTime, params.y = crtIntensity
    bg: vec4<f32>,
    fg: vec4<f32>,
    accent: vec4<f32>,
    response: vec4<f32>,
    params: vec4<f32>,
};

@group(0) @binding(0) var sceneSampler: sampler;
@group(0) @binding(1) var sceneTexture: texture_2d<f32>;
@group(0) @binding(2) var<uniform> theme: Theme;
@group(0) @binding(3) var prevFrameTexture: texture_2d<f32>;
@group(0) @binding(4) var orbitsTexture: texture_2d<f32>;


struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
};

@vertex
fn vertexMain(@builtin(vertex_index) vertexIndex: u32) -> VertexOutput {
    let pos = array<vec2<f32>, 6>(
        vec2<f32>(-1.0, -1.0), vec2<f32>( 1.0, -1.0), vec2<f32>(-1.0,  1.0),
        vec2<f32>(-1.0,  1.0), vec2<f32>( 1.0, -1.0), vec2<f32>( 1.0,  1.0)
    );
    var output: VertexOutput;
    let p = pos[vertexIndex];
    output.position = vec4<f32>(p, 0.0, 1.0);
    output.uv = (p * 0.5) + vec2<f32>(0.5, 0.5);
    output.uv.y = 1.0 - output.uv.y;
    return output;
}
// Simple present shader: sample the provided texture and output directly
@fragment
fn presentFragment(@location(0) uv: vec2<f32>) -> @location(0) vec4<f32> {
    let color = textureSample(sceneTexture, sceneSampler, uv).rgb;
    return vec4<f32>(color, 1.0);
}


@fragment
fn fragmentMain(@location(0) uv: vec2<f32>, @builtin(position) fragCoord: vec4<f32>) -> @location(0) vec4<f32> {
    // Barrel distortion warp
	var warpedUV = uv;
	let k1 = 0.1 * theme.params.y;
	let k2 = 0.02 * theme.params.y; // second-order term for physical lens model
	let center = vec2<f32>(0.5, 0.5);
	let radial = warpedUV - center;
	let r2 = dot(radial, radial);
	let distortion = 1.0 + k1 * r2 + k2 * (r2 * r2);
	warpedUV = center + radial * distortion;

    // Sample HDR scene color (raw HDR)
    let currentFrameHdr = textureSample(sceneTexture, sceneSampler, warpedUV).rgb;

    // Apply RGB response gain in HDR space.
    let rgb_energy = currentFrameHdr * theme.response.rgb;

    // Persistence operates in HDR: combine with decayed previous HDR frame.
    let prevFrameHdr = textureSample(prevFrameTexture, sceneSampler, warpedUV).rgb;
    let persistenceFactor = exp(-theme.params.x / 0.25);
    let persisted_hdr = rgb_energy + (prevFrameHdr * persistenceFactor);

    // Vignette in HDR before tonemapping
    var finalColorWithVignette = persisted_hdr;
    let vignetteStrength = 0.8 * theme.params.y;
    let vignetteFalloff = 0.5;
    let distFromCenter = length(uv - vec2<f32>(0.5, 0.5));
    let vignette = 1.0 - smoothstep(vignetteFalloff, 1.0, distFromCenter) * vignetteStrength;
    finalColorWithVignette *= vignette;

    // Sample orbits overlay with the same warp.
    let orbitsColor = textureSample(orbitsTexture, sceneSampler, warpedUV).rgb;
    finalColorWithVignette += orbitsColor;

    // Scanlines: dim every other row using integer row index
    let y = i32(floor(fragCoord.y));
    let scanlineFactor = 1.0 - (0.5 * theme.params.y);
    if ((y & 1) == 0) {
        finalColorWithVignette *= scanlineFactor;
    }

    // Tone map after persistence and vignette
    let tonemapped = finalColorWithVignette / (finalColorWithVignette + vec3<f32>(1.0));
    let ldrWithBg = tonemapped; // Background color addition disabled for now

    // Final safety: scrub NaN and clamp
    let noNan = select(vec3<f32>(0.0), ldrWithBg, ldrWithBg == ldrWithBg);
    let safeColor = clamp(noNan, vec3<f32>(0.0), vec3<f32>(1.0));
    return vec4<f32>(safeColor, 1.0);
}



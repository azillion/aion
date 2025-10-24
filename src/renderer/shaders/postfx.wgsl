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

struct OrbitMask {
    count: u32,
    targetRadiusPx: f32,
    screenSize: vec2<f32>,
    targetsPx: array<vec4<f32>, 16>,
};
@group(0) @binding(5) var<uniform> targetInfo: OrbitMask;

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

@fragment
fn fragmentMain(@location(0) uv: vec2<f32>, @builtin(position) fragCoord: vec4<f32>) -> @location(0) vec4<f32> {
    // Barrel distortion warp
    var warpedUV = uv;
    let barrelPower = 0.1 * theme.params.y;
    let center = vec2<f32>(0.5, 0.5);
    let radial = warpedUV - center;
    let r = length(radial);
    let distortion = 1.0 + barrelPower * r * r;
    warpedUV = center + radial / distortion;

    // Sample HDR scene color from the compute pass (raw HDR)
    let currentFrameHdr = textureSample(sceneTexture, sceneSampler, warpedUV).rgb;

    // Apply RGB response gain in HDR space.
    let rgb_energy = currentFrameHdr * theme.response.rgb;

    // Persistence operates in HDR: combine with decayed previous HDR frame.
    let prevFrameHdr = textureSample(prevFrameTexture, sceneSampler, warpedUV).rgb;
    let persistenceFactor = exp(-theme.params.x / 0.25);
    let persisted_hdr = rgb_energy + (prevFrameHdr * persistenceFactor);

    // Background glow (LDR contribution) will be added after tone mapping.
    var finalColorWithVignette = persisted_hdr;

    // Vignette
    let vignetteStrength = 0.8 * theme.params.y;
    let vignetteFalloff = 0.5;
    let distFromCenter = length(uv - vec2<f32>(0.5, 0.5));
    let vignette = 1.0 - smoothstep(vignetteFalloff, 1.0, distFromCenter) * vignetteStrength;
    finalColorWithVignette *= vignette;

    // Sample orbits overlay with the same warp.
    var orbitsColor = textureSample(orbitsTexture, sceneSampler, warpedUV).rgb;

    // Build a soft circular mask around bodies to make them readable over orbits
    var mask: f32 = 1.0;
    let fragPx = fragCoord.xy;
    let rPx = targetInfo.targetRadiusPx;
    // Smooth edge for better look
    let inner = rPx * 0.6;
    for (var i: u32 = 0u; i < targetInfo.count; i = i + 1u) {
        let c = targetInfo.targetsPx[i].xy;
        let d = distance(fragPx, c);
        let targetMask = smoothstep(inner, rPx, d);
        mask = min(mask, targetMask);
    }
    orbitsColor *= mask;
    finalColorWithVignette += orbitsColor;

    // Scanlines: dim every other row using integer row index
    let y = i32(floor(fragCoord.y));
    let scanlineFactor = 1.0 - (0.5 * theme.params.y);
    if ((y & 1) == 0) {
        finalColorWithVignette *= scanlineFactor;
    }

    // Tone map after persistence and vignette, then add subtle BG glow in LDR
    let tonemapped = finalColorWithVignette / (finalColorWithVignette + vec3<f32>(1.0));
    let ldrWithBg = tonemapped; // Temporarily disable background color addition

    // Final safety: scrub NaN and clamp
    let noNan = select(vec3<f32>(0.0), ldrWithBg, ldrWithBg == ldrWithBg);
    let safeColor = clamp(noNan, vec3<f32>(0.0), vec3<f32>(1.0));
    return vec4<f32>(safeColor, 1.0);
}

// Simple present shader: sample the provided texture and output directly
@fragment
fn presentFragment(@location(0) uv: vec2<f32>) -> @location(0) vec4<f32> {
    let color = textureSample(sceneTexture, sceneSampler, uv).rgb;
    return vec4<f32>(color, 1.0);
}



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

    // Sample HDR scene color
    let hdrColor = textureSample(sceneTexture, sceneSampler, warpedUV).rgb;

    // Luminance using theme response (smooth, no posterization)
    let luminance = dot(hdrColor, theme.response.rgb);

    // Phosphor persistence: sample previous frame and decay its luminance
    let prevFrameColor = textureSample(prevFrameTexture, sceneSampler, warpedUV).rgb;    
    // Guard against NaN/Inf in feedback
    let prevLen = length(prevFrameColor);
    let prevLuminance = clamp(prevLen / 1.732, 0.0, 1e6);
    let persistenceFactor = exp(-theme.params.x / 0.25);
    let persistedLuminance = max(luminance, prevLuminance * persistenceFactor);

    // Start with the base CRT background glow, scaled by intensity, then add scene phosphor
    let backgroundColor = theme.bg.rgb * theme.params.y;
    var finalColorWithVignette = backgroundColor + theme.fg.rgb * persistedLuminance;

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

    // Final safety: scrub NaN (x != x) and clamp huge values to avoid feedback loop
    let noNan = select(vec3<f32>(0.0), finalColorWithVignette, finalColorWithVignette == finalColorWithVignette);
    let safeColor = clamp(noNan, vec3<f32>(-1e6), vec3<f32>(1e6));
    return vec4<f32>(safeColor, 1.0);
}

// Simple present shader: sample the provided texture and output directly
@fragment
fn presentFragment(@location(0) uv: vec2<f32>) -> @location(0) vec4<f32> {
    let color = textureSample(sceneTexture, sceneSampler, uv).rgb;
    return vec4<f32>(color, 1.0);
}



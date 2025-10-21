struct Theme {
    bg: vec3<f32>,
    fg: vec3<f32>,
    accent: vec3<f32>,
    response: vec3<f32>,
    deltaTime: f32,
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
    let barrelPower = 0.1;
    let center = vec2<f32>(0.5, 0.5);
    let radial = warpedUV - center;
    let r = length(radial);
    let distortion = 1.0 + barrelPower * r * r;
    warpedUV = center + radial / distortion;

    // Sample HDR scene color
    let hdrColor = textureSample(sceneTexture, sceneSampler, warpedUV).rgb;

    // Luminance using theme response (smooth, no posterization)
    let luminance = dot(hdrColor, theme.response);

    // Phosphor persistence: sample previous frame and decay its luminance
    let prevFrameColor = textureSample(prevFrameTexture, sceneSampler, warpedUV).rgb;
    let prevLuminance = dot(prevFrameColor, vec3<f32>(0.2126, 0.7152, 0.0722));
    let persistenceFactor = exp(-theme.deltaTime / 0.25);
    let persistedLuminance = max(luminance, prevLuminance * persistenceFactor);

    // Theme application using persisted luminance
    var finalColor = theme.fg * persistedLuminance + theme.bg * (1.0 - persistedLuminance);

    // Vignette
    let vignetteStrength = 0.8;
    let vignetteFalloff = 0.5;
    let distFromCenter = length(uv - vec2<f32>(0.5, 0.5));
    let vignette = 1.0 - smoothstep(vignetteFalloff, 1.0, distFromCenter) * vignetteStrength;
    var finalColorWithVignette = finalColor * vignette;

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
    if ((y & 1) == 0) {
        finalColorWithVignette *= 0.85;
    }

    return vec4<f32>(finalColorWithVignette, 1.0);
}

// Simple present shader: sample the provided texture and output directly
@fragment
fn presentFragment(@location(0) uv: vec2<f32>) -> @location(0) vec4<f32> {
    let color = textureSample(sceneTexture, sceneSampler, uv).rgb;
    return vec4<f32>(color, 1.0);
}



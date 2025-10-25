// 3D Simplex Noise, adapted for WGSL.
// Source: https://github.com/stegu/webgl-noise

fn mod289_v3(x: vec3<f32>) -> vec3<f32> {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

fn mod289_v4(x: vec4<f32>) -> vec4<f32> {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

fn mod289_f(x: f32) -> f32 {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

fn permute_v3(x: vec3<f32>) -> vec3<f32> {
    return mod289_v3(((x*34.0)+1.0)*x);
}

fn permute_v4(x: vec4<f32>) -> vec4<f32> {
    return mod289_v4(((x*34.0)+1.0)*x);
}

fn snoise(v: vec3<f32>) -> f32 {
    let C = vec2<f32>(1.0/6.0, 1.0/3.0) ;
    let D = vec4<f32>(0.0, 0.5, 1.0, 2.0);

    // First corner
    var i  = floor(v + dot(v, C.yyy) );
    let x0 =   v - i + dot(i, C.xxx) ;

    // Other corners
    let g = step(x0.yzx, x0.xyz);
    let l = 1.0 - g;
    let i1 = min( g.xyz, l.zxy );
    let i2 = max( g.xyz, l.zxy );

    let x1 = x0 - i1 + C.xxx;
    let x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
    let x3 = x0 - D.yyy;      // -1.0+3.0*C.x = -0.5 = -D.y

    // Permutations
    i = mod289_v3(i);
    let p = permute_v4(
        permute_v4(vec4<f32>(i.y) + vec4<f32>(0.0, i1.y, i2.y, 1.0))
        + vec4<f32>(i.x) + vec4<f32>(0.0, i1.x, i2.x, 1.0)
    );

    // Gradients: 7x7 points over a rhombic dodecahedron.
    let n_ = 0.142857142857; // 1.0/7.0
    let ns = n_ * D.wyz - D.xzx;

    let j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)

    let x_ = floor(j * ns.z);
    let y_ = floor(j - 7.0 * x_ );    // mod(j,N)

    let x = x_ * ns.xxxx + ns.yyyy;
    let y = y_ * ns.xxxx + ns.yyyy;
    let h = 1.0 - abs(x) - abs(y);

    let b0 = vec4<f32>( x.xy, y.xy );
    let b1 = vec4<f32>( x.zw, y.zw );

    let s0 = floor(b0)*2.0 + 1.0;
    let s1 = floor(b1)*2.0 + 1.0;
    let sh = -step(h, vec4<f32>(0.0));

    let a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
    let a1 = b1.xzyw + s1.xzyw*sh.zzww ;

    var p0 = vec3<f32>(a0.xy,h.x);
    var p1 = vec3<f32>(a0.zw,h.y);
    var p2 = vec3<f32>(a1.xy,h.z);
    var p3 = vec3<f32>(a1.zw,h.w);

    // Normalise gradients
    let norm = inverseSqrt(vec4<f32>(dot(p0,p0), dot(p1,p1), dot(p2,p2), dot(p3,p3)));
    p0 = p0 * norm.x;
    p1 = p1 * norm.y;
    p2 = p2 * norm.z;
    p3 = p3 * norm.w;

    // Mix final noise value
    var m = max(0.6 - vec4<f32>(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), vec4<f32>(0.0));
    m = m * m;
    return 42.0 * dot( m*m, vec4<f32>( dot(p0,x0), dot(p1,x1), dot(p2,x2), dot(p3,x3) ) );
}

// --- Utility integer hash functions (PCG) and helpers ---
// Source references:
// - https://www.pcg-random.org/

fn pcg(n: u32) -> u32 {
    var h = n * 747796405u + 2891336453u;
    h = ((h >> ((h >> 28u) + 4u)) ^ h) * 277803737u;
    return (h >> 22u) ^ h;
}

fn pcg3d(p: vec3<u32>) -> vec3<u32> {
    var v = p * 1664525u + 1013904223u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    v ^= (v >> vec3<u32>(16u));
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    return v;
}

fn unorm_from_u32(x: u32) -> f32 { return f32(x) / 4294967295.0; }

fn rand01(n: u32) -> f32 { return unorm_from_u32(pcg(n)); }

fn rand3_from_u32(p: vec3<u32>) -> vec3<f32> {
    let h = pcg3d(p);
    return vec3<f32>(unorm_from_u32(h.x), unorm_from_u32(h.y), unorm_from_u32(h.z));
}


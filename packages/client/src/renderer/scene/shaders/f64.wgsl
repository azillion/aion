// Library for emulating 64-bit floating point arithmetic using two 32-bit floats.
// f64 Emulation Library

// Algorithm from "Accurate Sum and Dot Product" by T. Ogita, S. Rump, and S. Oishi
fn two_sum(a: f32, b: f32) -> vec2<f32> {
    let s = a + b;
    let t = (a - (s - b)) + (b - (s - a));
    return vec2<f32>(s, t);
}

fn two_sub(a: f32, b: f32) -> vec2<f32> {
    let s = a - b;
    let t = (a - (s + b)) + (b - (a - s));
    return vec2<f32>(s, t);
}

// Splits a single f32 into two f32s for Dekker multiplication
fn split(a: f32) -> vec2<f32> {
    let c = 134217729.0 * a; // 2^27 + 1
    let a_hi = c - (c - a);
    let a_lo = a - a_hi;
    return vec2<f32>(a_hi, a_lo);
}

// Emulated f64 multiplication
fn two_prod(a: f32, b: f32) -> vec2<f32> {
    let s = a * b;
    let a_split = split(a);
    let b_split = split(b);
    let err = (a_split.x * b_split.x - s) + (a_split.x * b_split.y) + (a_split.y * b_split.x) + (a_split.y * b_split.y);
    return vec2<f32>(s, err);
}

// Adds two emulated f64 numbers (a_h, a_l) + (b_h, b_l)
fn add_f64(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
    let s = two_sum(a.x, b.x);
    let e = two_sum(a.y, b.y);
    let c = s.y + e.x;
    let vh = s.x + c;
    let th = vh - s.x;
    let tl = c - th;
    let vl = tl + e.y;
    return vec2<f32>(vh, vl);
}

// Subtracts two emulated f64 numbers (a_h, a_l) - (b_h, b_l)
fn sub_f64(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
    let s = two_sub(a.x, b.x);
    let e = two_sub(a.y, b.y);
    let c = s.y + e.x;
    let vh = s.x + c;
    let th = vh - s.x;
    let tl = c - th;
    let vl = tl + e.y;
    return vec2<f32>(vh, vl);
}

// Subtracts two f64 (split into high/low pairs) and returns a single,
// high-precision f32 result. This is the primary function for calculating
// camera-relative positions.
fn sub_f64_to_f32(a_h: vec3<f32>, a_l: vec3<f32>, b_h: vec3<f32>, b_l: vec3<f32>) -> vec3<f32> {
    let r_x = sub_f64(vec2<f32>(a_h.x, a_l.x), vec2<f32>(b_h.x, b_l.x));
    let r_y = sub_f64(vec2<f32>(a_h.y, a_l.y), vec2<f32>(b_h.y, b_l.y));
    let r_z = sub_f64(vec2<f32>(a_h.z, a_l.z), vec2<f32>(b_h.z, b_l.z));
    return vec3<f32>(r_x.x + r_x.y, r_y.x + r_y.y, r_z.x + r_z.y);
}


// Adds an emulated f64 (a) and a standard f32 (b)
fn add_f64_f32(a: vec2<f32>, b: f32) -> vec2<f32> {
    // Fast path when high part is small enough that plain addition is safe
    if (abs(a.x) < 10000.0) {
        return vec2<f32>(a.x + b, a.y);
    }
    let s = two_sum(a.x, b);
    let vh = s.x;
    let vl = a.y + s.y;
    return vec2<f32>(vh, vl);
}



// Centralized Camera Uniform definition.
// This must match the layout expected by all shaders and populated by the CPU.
// Total size: 320 bytes.

struct CameraUniforms {
    // --- Matrices (64 bytes each) ---
    view: mat4x4<f32>,
    projection: mat4x4<f32>,
    viewProjection: mat4x4<f32>,

    // --- Vectors & Scalars (16 bytes per row) ---
    eye: vec3<f32>,
    _pad1: f32,

    forward: vec3<f32>,
    distance_to_target: f32,

    right: vec3<f32>,
    _pad2: f32,

    up: vec3<f32>,
    _pad3: f32,

    // NEW: Explicit vec4 for projection related constants
    projection_constants: vec4<f32>, // .x = lod_constant (height/vfov_rad)
};



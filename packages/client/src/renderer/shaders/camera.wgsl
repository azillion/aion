// Centralized Camera Uniform definition.
// This must match the layout populated by renderer/index.ts
// Total size: 352 bytes (88 floats).

struct CameraUniforms {
    // --- Matrices (64 bytes each) ---
    view: mat4x4<f32>,
    projection: mat4x4<f32>,
    viewProjection: mat4x4<f32>,

    // --- Vectors & Scalars (16 bytes per row) ---
    eye_pos_high: vec3<f32>,
    _pad_high: f32,

    eye_pos_low: vec3<f32>,
    _pad_low: f32,

    eye: vec3<f32>, // NOTE: This is (0,0,0) for camera-relative rendering
    _pad1: f32,

    forward: vec3<f32>,
    distance_to_look_at: f32,

    right: vec3<f32>,
    _pad2: f32, // This slot is unused in TS buffer

    up: vec3<f32>,
    _pad3: f32, // This slot is unused in TS buffer

    // .x = lod_constant, .y, .z, .w are unused for now
    projection_constants: vec4<f32>,
    
    // Final padding to reach 88 floats
    _pad4: vec4<f32>,
    _pad5: vec4<f32>,
    _pad6: vec4<f32>,
};



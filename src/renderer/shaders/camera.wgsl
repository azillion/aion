// Centralized Camera Uniform definition.
// This must match the layout expected by all shaders and populated by the CPU.
// Total size: 256 bytes.

struct CameraUniforms {
    // --- Matrices (64 bytes each) ---
    view: mat4x4<f32>,                // Offset: 0
    projection: mat4x4<f32>,          // Offset: 64
    viewProjection: mat4x4<f32>,      // Offset: 128

    // --- Vectors & Scalars (16 bytes per row) ---
    // We use vec3<f32> followed by f32 to tightly pack data into 16-byte alignment rows.
    
    // Row 1: Camera position
    eye: vec3<f32>,                   // Offset: 192
    _pad1: f32,                       // Offset: 204 (unused for now)

    // Row 2: Primary look direction (normalized)
    forward: vec3<f32>,               // Offset: 208
    distance_to_target: f32,          // Offset: 220 (Used for raytracing focus/DOF)

    // Row 3: Local Right axis (normalized)
    right: vec3<f32>,                 // Offset: 224
    _pad2: f32,                       // Offset: 236

    // Row 4: Local Up axis (normalized)
    up: vec3<f32>,                    // Offset: 240
    // Repurpose _pad3 to store viewport_height / vfov_rad for LOD calculations
    lod_constant: f32,                // Offset: 252
};



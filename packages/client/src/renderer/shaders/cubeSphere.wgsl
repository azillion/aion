// Mathematical helper functions for mapping between a 3D direction vector
// and a 2D UV coordinate on one of the 6 faces of a cube map.

// Input: A normalized 3D direction vector.
// Output: vec4f(u, v, faceIndex, 0.0), where u,v are in [0, 1] range.
fn directionToCubeUV(dir: vec3f) -> vec4f {
  let abs_dir = abs(dir);
  var face_index: f32;
  var uv: vec2f;

  if (abs_dir.x >= abs_dir.y && abs_dir.x >= abs_dir.z) { // X face
    face_index = select(1.0, 0.0, dir.x > 0.0);
    uv = vec2f(-dir.z, -dir.y) / abs_dir.x;
  } else if (abs_dir.y >= abs_dir.x && abs_dir.y >= abs_dir.z) { // Y face
    face_index = select(3.0, 2.0, dir.y > 0.0);
    uv = vec2f(dir.x, dir.z) / abs_dir.y;
  } else { // Z face
    face_index = select(5.0, 4.0, dir.z > 0.0);
    uv = vec2f(dir.x, -dir.y) / abs_dir.z;
  }
  return vec4f((uv * 0.5) + 0.5, face_index, 0.0);
}

// Input: uv in [0,1], faceIndex (0-5 integer).
// Output: A normalized 3D direction vector.
fn cubeUVToDirection(uv_in: vec2f, face_index: u32) -> vec3f {
    let uv = uv_in * 2.0 - 1.0; // Map uv from [0,1] to [-1,1]
    var dir: vec3f;
    switch (face_index) {
        case 0u: { dir = vec3f(1.0, -uv.y, -uv.x); } // +X
        case 1u: { dir = vec3f(-1.0, -uv.y, uv.x); }  // -X
        case 2u: { dir = vec3f(uv.x, 1.0, uv.y); }   // +Y
        case 3u: { dir = vec3f(uv.x, -1.0, -uv.y); } // -Y
        case 4u: { dir = vec3f(uv.x, -uv.y, 1.0); }   // +Z
        case 5u: { dir = vec3f(-uv.x, -uv.y, -1.0); } // -Z
        default: { dir = vec3f(1.0, 0.0, 0.0); }
    }
    return normalize(dir);
}

// Given a texel coordinate, face, cardinal direction offset, and texture width,
// return the neighbor's coordinate (in the same texel space) and face index.
fn getNeighbor(p: vec2i, face: u32, dir: vec2i, width: u32) -> vec3f {
    var np = vec2f(p) + vec2f(dir);
    var nf = f32(face);
    let W = f32(width);
    let max_coord = W - 1.0;

    if (np.x >= 0.0 && np.x < W && np.y >= 0.0 && np.y < W) {
        return vec3f(np, nf);
    }

    // Face indices: 0=+X, 1=-X, 2=+Y, 3=-Y, 4=+Z, 5=-Z
    let u = f32(p.x);
    let v = f32(p.y);

    if (np.y < 0.0) { // Top edge
        switch(face) {
            case 0u: { nf = 2.0; np = vec2f(max_coord, u); }
            case 1u: { nf = 2.0; np = vec2f(0.0, max_coord - u); }
            case 2u: { nf = 5.0; np = vec2f(max_coord - u, max_coord); }
            case 3u: { nf = 4.0; np = vec2f(u, 0.0); }
            case 4u: { nf = 2.0; np = vec2f(u, max_coord); }
            case 5u: { nf = 2.0; np = vec2f(max_coord, max_coord - u); }
            default: {}
        }
    } else if (np.y >= W) { // Bottom edge
        switch(face) {
            case 0u: { nf = 3.0; np = vec2f(0.0, u); }
            case 1u: { nf = 3.0; np = vec2f(max_coord, max_coord - u); }
            case 2u: { nf = 4.0; np = vec2f(u, 0.0); }
            case 3u: { nf = 5.0; np = vec2f(max_coord - u, 0.0); }
            case 4u: { nf = 3.0; np = vec2f(u, max_coord); }
            case 5u: { nf = 3.0; np = vec2f(0.0, max_coord - u); }
            default: {}
        }
    } else if (np.x < 0.0) { // Left edge
        switch(face) {
            case 0u: { nf = 4.0; np = vec2f(v, 0.0); }
            case 1u: { nf = 5.0; np = vec2f(max_coord - v, 0.0); }
            case 2u: { nf = 1.0; np = vec2f(v, 0.0); }
            case 3u: { nf = 1.0; np = vec2f(max_coord - v, max_coord); }
            case 4u: { nf = 1.0; np = vec2f(max_coord, v); }
            case 5u: { nf = 1.0; np = vec2f(0.0, v); }
            default: {}
        }
    } else if (np.x >= W) { // Right edge
        switch(face) {
            case 0u: { nf = 5.0; np = vec2f(max_coord - v, max_coord); }
            case 1u: { nf = 4.0; np = vec2f(v, max_coord); }
            case 2u: { nf = 0.0; np = vec2f(max_coord - v, 0.0); }
            case 3u: { nf = 0.0; np = vec2f(v, max_coord); }
            case 4u: { nf = 0.0; np = vec2f(0.0, v); }
            case 5u: { nf = 0.0; np = vec2f(max_coord, v); }
            default: {}
        }
    }

    return vec3f(np, nf);
}



struct GridVertices {
  data: array<vec3f>,
};

struct GridElevations {
  data: array<f32>,
};

struct GridIndices {
  data: array<u32>,
};

@group(0) @binding(7) var<storage, read> grid_vertices: GridVertices;
@group(0) @binding(8) var<storage, read> grid_elevations: GridElevations;
@group(0) @binding(9) var<storage, read> grid_indices: GridIndices;

fn barycentric_weights(p: vec3f, a: vec3f, b: vec3f, c: vec3f) -> vec3f {
  let v0 = b - a;
  let v1 = c - a;
  let v2 = p - a;
  let d00 = dot(v0, v0);
  let d01 = dot(v0, v1);
  let d11 = dot(v1, v1);
  let d20 = dot(v2, v0);
  let d21 = dot(v2, v1);
  let denom = d00 * d11 - d01 * d01;
  let v = (d11 * d20 - d01 * d21) / denom;
  let w = (d00 * d21 - d01 * d20) / denom;
  let u = 1.0 - v - w;
  return vec3f(u, v, w);
}

fn project_point_to_triangle_plane(p: vec3f, a: vec3f, b: vec3f, c: vec3f) -> vec3f {
  let n = normalize(cross(b - a, c - a));
  let d = dot(n, p - a);
  return p - d * n;
}

fn inside_barycentric(bc: vec3f) -> bool {
  return bc.x >= 0.0 && bc.y >= 0.0 && bc.z >= 0.0 && (abs(bc.x + bc.y + bc.z - 1.0) <= 1e-3);
}

// Returns interpolated elevation sampled from the coarse grid
fn interpolate_triangle(p: vec3f, a: vec3f, b: vec3f, c: vec3f, ea: f32, eb: f32, ec: f32) -> f32 {
  let proj = project_point_to_triangle_plane(p, a, b, c);
  let bc = barycentric_weights(proj, a, b, c);
  if (!inside_barycentric(bc)) { return -1.0; }
  return bc.x * ea + bc.y * eb + bc.z * ec;
}

// Slow reference implementation: linear scan over all triangles
fn sampleCoarseGrid(world_pos: vec3f) -> f32 {
  let tri_count = arrayLength(&grid_indices.data) / 3u;
  for (var i: u32 = 0u; i < tri_count; i = i + 1u) {
    let ia = grid_indices.data[i * 3u + 0u];
    let ib = grid_indices.data[i * 3u + 1u];
    let ic = grid_indices.data[i * 3u + 2u];
    let a = grid_vertices.data[ia];
    let b = grid_vertices.data[ib];
    let c = grid_vertices.data[ic];
    let ea = grid_elevations.data[ia];
    let eb = grid_elevations.data[ib];
    let ec = grid_elevations.data[ic];
    let h = interpolate_triangle(world_pos, a, b, c, ea, eb, ec);
    if (h >= 0.0) { return h; }
  }
  return 0.0;
}



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

// NEW Barycentric weights for a point p inside a spherical triangle (a,b,c)
// All vectors are assumed to be normalized (directions from the sphere center).
fn barycentric_weights_spherical(p: vec3f, a: vec3f, b: vec3f, c: vec3f) -> vec3f {
  let n = cross(b - a, c - a); // Normal of the planar triangle
  // Areas of the sub-triangles
  let area_pbc = dot(n, cross(b - p, c - p));
  let area_pca = dot(n, cross(c - p, a - p));
  let area_pab = dot(n, cross(a - p, b - p));
  // The total area is the sum of the sub-areas.
  let total_area_inv = 1.0 / (area_pbc + area_pca + area_pab);
  // The weights are the ratio of the sub-areas to the total area.
  return vec3f(area_pbc, area_pca, area_pab) * total_area_inv;
}

// NEW interpolation function using spherical logic.
fn interpolate_triangle(p: vec3f, a: vec3f, b: vec3f, c: vec3f, ea: f32, eb: f32, ec: f32) -> f32 {
  // Determine the triangle's winding order. The cross product of two edges
  // gives us a vector that points "outward" from the triangle's plane.
  let outward_normal = normalize(cross(b - a, c - a));

  // A point is inside the spherical triangle if it is "behind" all three
  // of the great-circle planes defined by the edges. We can test this by
  // checking the sign of the dot product between the point and the plane's normal.
  // The plane normal for a great-circle arc (e.g., from A to B) is simply cross(A, B).

  // We determine the correct "sign" by checking which way the outward_normal points
  // relative to the center of the sphere. The vectors a,b,c point from the center.
  let winding = sign(dot(outward_normal, a));

  if (sign(dot(cross(a, b), p)) != winding) { return -1.0; }
  if (sign(dot(cross(b, c), p)) != winding) { return -1.0; }
  if (sign(dot(cross(c, a), p)) != winding) { return -1.0; }

  // If all checks pass, the point is unambiguously inside.
  let bc = barycentric_weights_spherical(p, a, b, c);
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
    if (h >= 0.0) { 
        return h; // Return just the height
    }
  }
  return 0.0; // Fallback
}

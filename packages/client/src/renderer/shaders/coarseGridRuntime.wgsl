@group(0) @binding(11) var terrain_height_tex: texture_2d_array<f32>;

fn sampleCoarseGrid(world_dir: vec3f) -> f32 {
  let dir = normalize(world_dir);
  let cube = directionToCubeUV(dir);
  let dims = textureDimensions(terrain_height_tex);
  let tex_size = vec2f(f32(dims.x), f32(dims.y));
  if (tex_size.x < 1.0 || tex_size.y < 1.0) {
    return 0.0;
  }

  let face = i32(clamp(round(cube.z), 0.0, 5.0));
  let coord = cube.xy * (tex_size - vec2f(1.0));
  let base = floor(coord);
  let frac = coord - base;

  let base_i_f = clamp(base, vec2f(0.0), tex_size - vec2f(1.0));
  let base_j_f = clamp(base + vec2f(1.0), vec2f(0.0), tex_size - vec2f(1.0));
  let base_i = vec2i(i32(base_i_f.x), i32(base_i_f.y));
  let base_j = vec2i(i32(base_j_f.x), i32(base_j_f.y));

  let h00 = textureLoad(terrain_height_tex, base_i, face, 0).r;
  let h10 = textureLoad(terrain_height_tex, vec2i(base_j.x, base_i.y), face, 0).r;
  let h01 = textureLoad(terrain_height_tex, vec2i(base_i.x, base_j.y), face, 0).r;
  let h11 = textureLoad(terrain_height_tex, base_j, face, 0).r;

  let hx0 = mix(h00, h10, frac.x);
  let hx1 = mix(h01, h11, frac.x);
  return mix(hx0, hx1, frac.y);
}



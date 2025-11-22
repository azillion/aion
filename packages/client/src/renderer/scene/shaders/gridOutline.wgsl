#include "camera.wgsl"

struct PlanetUniform {
  position: vec3f,
  radius: f32,
};

@group(0) @binding(0) var<uniform> camera: CameraUniforms;
@group(0) @binding(1) var<uniform> planet: PlanetUniform;

struct VSOut {
  @builtin(position) clip_pos: vec4f,
};

@vertex
fn vs_main(@location(0) pos: vec3f) -> VSOut {
  let world = planet.position + pos * (planet.radius * 1.001);
  var out: VSOut;
  out.clip_pos = camera.viewProjection * vec4f(world, 1.0);
  return out;
}

@fragment
fn fs_main() -> @location(0) vec4f {
  return vec4f(0.0, 0.0, 0.0, 1.0);
}



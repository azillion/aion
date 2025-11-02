import type { Body } from '@shared/types';

const FLOATS_PER_BODY_GPU = 32;

export interface ProcessedSceneData {
  data: Float32Array;
  count: number;
}

/**
 * Serializes tiered Body arrays into zero-padded Float32Arrays matching
 * the WGSL Sphere struct layout. Eliminates per-frame object churn.
 */
export class SceneDataProcessor {
  /**
   * Splits a 64-bit float (a JS `number`) into two 32-bit floats (high and low)
   * for emulated double-precision arithmetic on the GPU.
   * @param v The 64-bit float to split.
   * @returns A tuple `[high, low]` where `high` is the high-order bits and `low` is the remainder.
   */
  public splitDouble(v: number): [number, number] {
    if (v === 0.0) return [0, 0];
    const high = Math.fround(v);
    const low = v - high;
    return [high, low];
  }

  public process(bodies: Body[]): ProcessedSceneData {
    const buffer = new Float32Array(bodies.length * FLOATS_PER_BODY_GPU);
    this.serializeBodies(buffer, bodies);
    return { data: buffer, count: bodies.length };
  }

  private serializeBodies(buffer: Float32Array, bodies: Body[]) {
    for (let i = 0; i < bodies.length; i++) {
      const body = bodies[i];
      const base = i * FLOATS_PER_BODY_GPU;

      // Split the world-space f64 position into high and low f32 parts.
      const [posX_h, posX_l] = this.splitDouble(body.position[0]);
      const [posY_h, posY_l] = this.splitDouble(body.position[1]);
      const [posZ_h, posZ_l] = this.splitDouble(body.position[2]);

      // pos_high_and_radius: vec4<f32>
      buffer[base + 0] = posX_h;
      buffer[base + 1] = posY_h;
      buffer[base + 2] = posZ_h;
      buffer[base + 3] = body.radius;

      // pos_low_and_unused: vec4<f32>
      buffer[base + 4] = posX_l;
      buffer[base + 5] = posY_l;
      buffer[base + 6] = posZ_l;
      // buffer[base + 7] is unused

      // albedo_and_atmos_flag: vec4<f32>
      buffer.set(body.albedo, base + 8);
      buffer[base + 11] = body.terrain?.atmosphere ? 1.0 : 0.0;

      // emissive_and_terrain_flag: vec4<f32>
      const emissive = body.emissive ?? [0, 0, 0];
      buffer.set(emissive, base + 12);
      buffer[base + 15] = body.terrain ? 1.0 : 0.0;

      // ref_idx_opacity_pad: vec4<f32>
      buffer[base + 16] = 1.0; // ref_idx
      buffer[base + 17] = 1.0; // opacity
      // buffer[base + 18, 19] are unused

      // terrain_params: vec4<f32>
      if (body.terrain) {
        buffer[base + 20] = body.terrain.radius;
        buffer[base + 21] = body.terrain.seaLevel;
        buffer[base + 22] = body.terrain.maxHeight;
        buffer[base + 23] = body.terrain.noiseSeed;
      } else {
        buffer.fill(0.0, base + 20, base + 24);
      }
      // The rest of the buffer (padding1, padding2) is left as zeros.
    }
  }
}



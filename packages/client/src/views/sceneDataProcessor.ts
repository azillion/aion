import type { Body } from '@shared/types';

const FLOATS_PER_BODY_GPU = 24; // Must match the Sphere struct in WGSL

export interface ProcessedSceneData {
  nearData: Float32Array;
  midData: Float32Array;
  farData: Float32Array;
  nearCount: number;
  midCount: number;
  farCount: number;
}

/**
 * Serializes tiered Body arrays into zero-padded Float32Arrays matching
 * the WGSL Sphere struct layout. Eliminates per-frame object churn.
 */
export class SceneDataProcessor {
  public process(tieredScene: { near: Body[]; mid: Body[]; far: Body[]; }): ProcessedSceneData {
    // Fresh, perfectly-sized buffers each call for robustness
    const nearBuffer = new Float32Array(tieredScene.near.length * FLOATS_PER_BODY_GPU);
    const midBuffer = new Float32Array(tieredScene.mid.length * FLOATS_PER_BODY_GPU);
    const farBuffer = new Float32Array(tieredScene.far.length * FLOATS_PER_BODY_GPU);

    this.serializeTier(nearBuffer, tieredScene.near);
    this.serializeTier(midBuffer, tieredScene.mid);
    this.serializeTier(farBuffer, tieredScene.far);

    return {
      nearData: nearBuffer,
      midData: midBuffer,
      farData: farBuffer,
      nearCount: tieredScene.near.length,
      midCount: tieredScene.mid.length,
      farCount: tieredScene.far.length,
    };
  }

  private serializeTier(buffer: Float32Array, bodies: Body[]) {
    for (let i = 0; i < bodies.length; i++) {
      const body = bodies[i];
      const base = i * FLOATS_PER_BODY_GPU;

      buffer[base + 0] = body.position[0];
      buffer[base + 1] = body.position[1];
      buffer[base + 2] = body.position[2];
      buffer[base + 3] = body.radius;

      buffer.set(body.albedo, base + 4);
      buffer[base + 7] = body.terrain?.atmosphere ? 1.0 : 0.0;

      const emissive = body.emissive ?? [0, 0, 0];
      buffer.set(emissive, base + 8);
      buffer[base + 11] = body.terrain ? 1.0 : 0.0;

      buffer[base + 12] = 1.0; // ref_idx
      buffer[base + 13] = 1.0; // opacity
      buffer[base + 14] = 0.0;
      buffer[base + 15] = 0.0;

      if (body.terrain) {
        buffer[base + 16] = body.terrain.radius;
        buffer[base + 17] = body.terrain.seaLevel;
        buffer[base + 18] = body.terrain.maxHeight;
        buffer[base + 19] = body.terrain.noiseSeed;
      } else {
        buffer.fill(0.0, base + 16, base + 20);
      }
      // base + 20..23 padding left as zeros
    }
  }
}



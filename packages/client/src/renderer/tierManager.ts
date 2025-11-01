import type { Body, SystemState, Vec3 } from '@shared/types';
import { NEAR_TIER_CUTOFF, MID_TIER_CUTOFF, MID_TIER_SCALE, FAR_TIER_SCALE } from '@shared/constants';

export interface TieredScene {
  near: Body[];
  mid: Body[];
  far: Body[];
  shadowCasters: Body[];
}

export class TierManager {
  public build(systemState: SystemState, cameraWorldPosition: Vec3, playerShipId: string | null): TieredScene {
    const near: Body[] = [];
    const mid: Body[] = [];
    const far: Body[] = [];
    const shadowCasters: Body[] = [];

    systemState.bodies.forEach(body => {
      if (body.id === playerShipId) return;

      const dx = body.position[0] - cameraWorldPosition[0];
      const dy = body.position[1] - cameraWorldPosition[1];
      const dz = body.position[2] - cameraWorldPosition[2];
      const dist_to_center = Math.sqrt(dx*dx + dy*dy + dz*dz);
      const dist_to_surface = Math.max(0, dist_to_center - body.radius);

      if (body.mass > 1e22) {
        const hasAtmosphere = body.terrain?.atmosphere !== undefined;
        const ATMOSPHERE_RADIUS_SCALE = 1.025;
        let shadowRadius = body.terrain ? body.terrain.radius : body.radius;
        if (hasAtmosphere) { shadowRadius *= ATMOSPHERE_RADIUS_SCALE; }
        shadowCasters.push({ ...body, position: [dx, dy, dz], radius: shadowRadius });
      }

      if (dist_to_surface < NEAR_TIER_CUTOFF) {
        near.push({ ...body, position: [dx, dy, dz] });
      } else if (dist_to_surface < MID_TIER_CUTOFF) {
        mid.push({ ...body, position: [dx / MID_TIER_SCALE, dy / MID_TIER_SCALE, dz / MID_TIER_SCALE], radius: body.radius / MID_TIER_SCALE });
      } else {
        far.push({ ...body, position: [dx / FAR_TIER_SCALE, dy / FAR_TIER_SCALE, dz / FAR_TIER_SCALE], radius: body.radius / FAR_TIER_SCALE });
      }
    });

    return { near, mid, far, shadowCasters };
  }
}



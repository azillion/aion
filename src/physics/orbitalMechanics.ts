import { vec3 } from 'gl-matrix';
import { G } from '../shared/constants';
import type { Vec3 } from '../shared/types';

const ORBIT_SAMPLES = 256;
export const ORBIT_MAX_POINTS = ORBIT_SAMPLES + 1;

export interface OrbitalElements {
  a: number; // semi-major axis
  e: number; // eccentricity
  currentTrueAnomaly: number;
  points: Float32Array;
  pointCount: number;
  periapsis: Vec3 | null;
  apoapsis: Vec3 | null;
  ascendingNode: Vec3 | null;
  descendingNode: Vec3 | null;
}

export function calculateAnalyticOrbit(
  r_vec: Vec3,
  v_vec: Vec3,
  totalSystemMass: number,
  systemScale: number
): OrbitalElements | null {
    const r = Math.hypot(...r_vec);
    const v = Math.hypot(...v_vec);
    const mu = G * totalSystemMass;
    const h_vec = vec3.cross(vec3.create(), r_vec, v_vec) as Vec3;
    const h = vec3.length(h_vec);
    if (h < 1e-9) return null; // No angular momentum, not an orbit

    const c_vec = vec3.cross(vec3.create(), v_vec, h_vec) as Vec3;
    const r_norm = vec3.normalize(vec3.create(), r_vec);
    
    const e_vec_term1 = vec3.scale(vec3.create(), c_vec, 1 / mu);
    const e_vec = vec3.subtract(vec3.create(), e_vec_term1, r_norm) as Vec3;
    const e = vec3.length(e_vec);

    const invA = 2 / r - (v * v) / mu;
    if (!isFinite(invA) || invA <= 0 || e >= 1) {
      return null;
    }
    const a = 1 / invA;
    const bSemi = a * Math.sqrt(Math.max(0, 1 - e * e));

    const W = vec3.normalize(vec3.create(), h_vec) as Vec3;
    let P_unnormalized = e > 1e-6 ? e_vec : r_vec;
    let P = vec3.normalize(vec3.create(), P_unnormalized) as Vec3;

    // Project P onto the orbital plane to ensure orthogonality with W
    const dotPW = vec3.dot(P, W);
    P = vec3.scaleAndAdd(vec3.create(), P, W, -dotPW) as Vec3;
    vec3.normalize(P, P);

    const Q = vec3.cross(vec3.create(), W, P) as Vec3;

    // Generate orbit ribbon points
    const points = new Float32Array(ORBIT_MAX_POINTS * 4);
    for (let k = 0; k < ORBIT_MAX_POINTS; k++) {
      const theta = (k / (ORBIT_MAX_POINTS - 1)) * Math.PI * 2;
      const xPeri = a * (Math.cos(theta) - e);
      const yPeri = bSemi * Math.sin(theta);
      const base = k * 4;

      const pos = vec3.create();
      vec3.scaleAndAdd(pos, pos, P, xPeri);
      vec3.scaleAndAdd(pos, pos, Q, yPeri);
      vec3.scale(pos, pos, systemScale);
      
      points[base+0] = pos[0];
      points[base+1] = pos[1];
      points[base+2] = pos[2];
      points[base+3] = 1.0;
    }
    
    const dotR_Q = vec3.dot(r_vec, Q);
    const dotR_P = vec3.dot(r_vec, P);
    const currentTrueAnomaly = Math.atan2(dotR_Q, dotR_P);

    let periapsis: Vec3 | null = null;
    let apoapsis: Vec3 | null = null;
    if (e > 1e-6) {
        periapsis = vec3.scale(vec3.create(), P, a * (1.0 - e) * systemScale) as Vec3;
        apoapsis = vec3.scale(vec3.create(), P, -a * (1.0 + e) * systemScale) as Vec3;
    }

    let ascendingNode: Vec3 | null = null;
    let descendingNode: Vec3 | null = null;
    const nodeLine: Vec3 = [-h_vec[1], h_vec[0], 0];
    if (vec3.length(nodeLine) > 1e-9) {
        const n = vec3.normalize(vec3.create(), nodeLine);
        const nu = Math.atan2(vec3.dot(n, Q), vec3.dot(n, P));
        const r_dist = a * (1 - e * e) / (1 + e * Math.cos(nu));
        ascendingNode = vec3.scale(vec3.create(), n, r_dist * systemScale) as Vec3;
        descendingNode = vec3.scale(vec3.create(), n, -r_dist * systemScale) as Vec3;
    }

    return {
        a, e, currentTrueAnomaly,
        points, pointCount: ORBIT_MAX_POINTS,
        periapsis, apoapsis, ascendingNode, descendingNode,
    };
}



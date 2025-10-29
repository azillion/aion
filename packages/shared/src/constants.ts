// Tier constants for multi-scale rendering
export const NEAR_TIER_CUTOFF = 1000.0;     // 1,000 km. The high-precision zone for terrain and local objects.
export const MID_TIER_CUTOFF = 1.0e6;       // 1 million km. The local space zone.

// Scale factors to convert world units (km) to local tier units
// By scaling coordinates down, we keep them in a high-precision range for f32 numbers.
export const MID_TIER_SCALE = 1000.0;       // 1 mid-tier unit = 1,000 km
export const FAR_TIER_SCALE = 100000.0;     // 1 far-tier unit = 100,000 km

/**
 * The universal gravitational constant.
 * Units: km^3 / (kg * s^2)
 */
export const G = 6.67430e-20;

// Planetary sphere-of-influence heuristic for switching to planet-centric reference frames
export const PLANETARY_SOI_RADIUS_MULTIPLIER = 100.0;



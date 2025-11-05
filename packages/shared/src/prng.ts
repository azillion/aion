/**
 * A simple and fast 32-bit pseudo-random number generator (PRNG)
 * using the Mulberry32 algorithm.
 */
export class PRNG {
  private seed: number;

  constructor(seed: number) {
    this.seed = seed;
  }

  /**
   * Returns a random float between 0 (inclusive) and 1 (exclusive).
   * This implementation uses the Mulberry32 algorithm, which has better
   * statistical properties than simple linear congruential generators.
   */
  public nextFloat(): number {
    let t = (this.seed += 0x6d2b79f5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  /**
   * Returns a random float within a given range [min, max).
   * @param min The minimum value (inclusive).
   * @param max The maximum value (exclusive).
   */
  public nextInRange(min: number, max: number): number {
    return this.nextFloat() * (max - min) + min;
  }
}
export class PRNG {
  private seed: number;

  constructor(seed: number) {
    this.seed = seed;
  }

  // Returns a random float between 0 (inclusive) and 1 (exclusive)
  public nextFloat(): number {
    this.seed |= 0;
    this.seed = (this.seed + 0x9e3779b9) | 0;
    let t = this.seed ^ (this.seed >>> 16);
    t = Math.imul(t, 0x21f0aaad);
    t = t ^ (t >>> 15);
    t = Math.imul(t, 0x735a2d97);
    return ((t = t ^ (t >>> 15)) >>> 0) / 4294967296;
  }
  
  // Returns a random float within a given range
  public nextInRange(min: number, max: number): number {
    return this.nextFloat() * (max - min) + min;
  }
}



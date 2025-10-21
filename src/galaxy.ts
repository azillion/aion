import { PRNG } from './shared/prng';
import type { Vec3, Star } from './shared/types';

export class Galaxy {
  public stars: Star[] = [];

  constructor(seed: number, numStars: number, cubeSize: number) {
    const prng = new PRNG(seed);
    const halfSize = cubeSize / 2;

    for (let i = 0; i < numStars; i++) {
      const position: Vec3 = [
        prng.nextInRange(-halfSize, halfSize),
        prng.nextInRange(-halfSize, halfSize),
        prng.nextInRange(-halfSize, halfSize)
      ];
      
      const colorR = 1.0 - prng.nextFloat() * 0.2;
      const colorG = 1.0 - prng.nextFloat() * 0.2;
      const colorB = 1.0;
      const color: Vec3 = [colorR, colorG, colorB];

      const size = prng.nextInRange(0.5, 2.0);

      this.stars.push({ position, color, size });
    }
  }
}

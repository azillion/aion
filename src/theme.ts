import type { Vec3 } from './shared/types';

export interface Theme {
  bg: Vec3;
  fg: Vec3;
  accent: Vec3;
}

export const themes: { [key: string]: Theme } = {
  amber: {
    bg: [0.165, 0.118, 0.063],
    fg: [0.961, 0.725, 0.420],
    accent: [0.988, 0.863, 0.643],
  },
  green: {
    bg: [0.031, 0.071, 0.031],
    fg: [0.494, 0.988, 0.522],
    accent: [0.694, 1.0, 0.745],
  },
  white: {
    bg: [0.039, 0.039, 0.039],
    fg: [0.910, 0.910, 0.910],
    accent: [1.0, 1.0, 1.0],
  },
};



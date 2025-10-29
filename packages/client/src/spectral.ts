import type { Vec3 } from '@shared/types';

export const spectralResponses: { [key: string]: Vec3 } = {
  'Full Color': [1.0, 1.0, 1.0],
  'Visible (Y)': [0.2126, 0.7152, 0.0722],
  'Infrared (IR)': [0.8, 0.15, 0.05],
  'Ultraviolet (UV)': [0.1, 0.2, 0.7],
  'Science (R-Band)': [1.0, 0.0, 0.0],
};



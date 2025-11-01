import type { Body, Vec3 } from '@shared/types';
import { PLANETARY_SOI_RADIUS_MULTIPLIER } from '@shared/constants';

export class ReferenceFrameManager {
  /**
   * Determines the most appropriate reference body based on the camera's proximity.
   * @param bodies - The list of all bodies in the system.
   * @param worldCameraEye - The camera's true position in world space.
   * @param currentReferenceBodyId - The ID of the currently active reference body.
   * @returns The ID of the new reference body if a switch should occur.
   */
  public update(
    bodies: Body[],
    worldCameraEye: Vec3,
    currentReferenceBodyId: string | null
  ): string | null {
    const currentReferenceBody = bodies.find(b => b.id === currentReferenceBodyId) ?? bodies.find(b => b.id === 'sol');
    let newReferenceBody = currentReferenceBody;
    let closestDistSq = Infinity;

    bodies.forEach(body => {
      if (body.mass < 1e22) return; // Only consider sufficiently massive bodies
      const dx = body.position[0] - worldCameraEye[0];
      const dy = body.position[1] - worldCameraEye[1];
      const dz = body.position[2] - worldCameraEye[2];
      const distSq = dx*dx + dy*dy + dz*dz;
      const soiRadius = body.radius * PLANETARY_SOI_RADIUS_MULTIPLIER;
      if (distSq < (soiRadius * soiRadius) && distSq < closestDistSq) {
        closestDistSq = distSq;
        newReferenceBody = body;
      }
    });

    if (closestDistSq === Infinity && currentReferenceBodyId !== 'sol') {
      newReferenceBody = bodies.find(b => b.id === 'sol');
    }

    if (newReferenceBody && currentReferenceBodyId !== newReferenceBody.id) {
      console.log(`Switching reference frame to: ${newReferenceBody.name}`);
      return newReferenceBody.id;
    }

    return null; // No change
  }
}



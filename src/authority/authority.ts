import type { SystemState } from '../shared/types';

export interface Authority {
	/**
	 * Queries the current state of the simulated system.
	 */
	query(): Promise<SystemState>;

	/**
	 * Advances the simulation by a given time step.
	 * @param deltaTime The time elapsed in seconds.
	 */
	tick(deltaTime: number): Promise<void>;
}

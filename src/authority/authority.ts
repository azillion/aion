import type { Body, SystemState } from '../shared/types';

export interface InputState { deltaX: number; deltaY: number; keys: Set<string>; }

export interface Authority {
	/**
	 * Queries the current state of the simulated system.
	 */
	query(): Promise<SystemState>;

	/**
	 * Advances the simulation by a given time step.
	 * @param deltaTime The time elapsed in seconds.
	 */
    tick(deltaTime: number, input: InputState): Promise<void>;

	/**
	 * Sets the simulation's time scale.
	 * @param scale A multiplier for the passage of time.
	 */
	setTimeScale(scale: number): void;

	/**
	 * Adds a new body to the simulation.
	 */
	addBody(body: Omit<Body, 'id'>): void;

	/**
	 * Initiates an auto-landing sequence for the player ship onto the target body.
	 * If targetBodyId is null or invalid, the call is ignored.
	 */
	autoLand(targetBodyId: string | null): void;

    /**
     * Instantly moves the player ship to just above the surface of the target body
     * and orients the ship so that its up vector matches the local horizon normal.
     * If targetBodyId is null or invalid, the call is ignored.
     */
    teleportToSurface(targetBodyId: string | null): void;
}

import type { Body, SystemState } from '@shared/types';

export interface InputState { keys: Set<string>; }

export interface Authority {
	query(): Promise<SystemState>;
	tick(deltaTime: number, input: InputState): Promise<void>;
	setTimeScale(scale: number): void;
	addBody(body: Omit<Body, 'id'>): void;
	autoLand(targetBodyId: string | null): void;
	teleportToSurface(targetBodyId: string | null): void;
}



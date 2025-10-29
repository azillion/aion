import type { Body, SystemState } from './types';

export interface InputState { keys: Set<string>; }

export interface Authority {
    query(): Promise<SystemState>;
    tick(deltaTime: number, input: InputState): Promise<void>;
    setTimeScale(scale: number): Promise<void>;
    addBody(body: Omit<Body, 'id'>): Promise<void>;
    autoLand(targetBodyId: string | null): Promise<void>;
    teleportToSurface(targetBodyId: string | null): Promise<void>;
}



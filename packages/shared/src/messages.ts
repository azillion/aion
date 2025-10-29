import type { InputState } from './authority';
import type { Body, SystemState } from './types';

// --- Client -> Server Messages ---

export interface QueryMessage {
  type: 'query';
  queryId: number;
}

export interface TickMessage {
  type: 'tick';
  deltaTime: number;
  input: InputState;
}

export interface SetTimeScaleMessage {
  type: 'setTimeScale';
  scale: number;
}

export interface AddBodyMessage {
  type: 'addBody';
  body: Omit<Body, 'id'>;
}

export interface AutoLandMessage {
  type: 'autoLand';
  targetBodyId: string | null;
}

export interface TeleportToSurfaceMessage {
  type: 'teleportToSurface';
  targetBodyId: string | null;
}

/** A union of all possible messages sent FROM the client TO the server. */
export type ClientToServerMessage =
  | QueryMessage
  | TickMessage
  | SetTimeScaleMessage
  | AddBodyMessage
  | AutoLandMessage
  | TeleportToSurfaceMessage;

// --- Server -> Client Messages ---

export interface QueryResultMessage {
  type: 'queryResult';
  queryId: number;
  state: SystemState;
}

/** A union of all possible messages sent FROM the server TO the client. */
export type ServerToClientMessage = QueryResultMessage;



import type { ClientToServerMessage, ServerToClientMessage } from '@shared/messages';

export interface IAuthorityConnection {
  postMessage(message: ClientToServerMessage): Promise<void>;
  onMessage: ((message: ServerToClientMessage) => void) | null;
}

export interface IAuthorityProvider {
  createConnection(): Promise<IAuthorityConnection>;
}



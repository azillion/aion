import type { Authority, InputState } from '@shared/authority';
import type { Body, SystemState } from '@shared/types';
import type { ClientToServerMessage, ServerToClientMessage } from '@shared/messages';
import type { QueryMessage } from '@shared/messages';

// This interface defines the contract for any transport layer (in-memory, websocket, etc.)
export interface IAuthorityConnection {
    postMessage(message: ClientToServerMessage): Promise<void>;
    onMessage: ((message: ServerToClientMessage) => void) | null;
}

export class ClientAuthority implements Authority {
    private connection: IAuthorityConnection;
    private pendingQueries: Map<number, (value: SystemState) => void> = new Map();
    private queryIdCounter = 0;

    constructor(connection: IAuthorityConnection) {
        this.connection = connection;
        this.connection.onMessage = this.handleMessage.bind(this);
    }

    private handleMessage(message: ServerToClientMessage): void {
        if (message.type === 'queryResult' && this.pendingQueries.has(message.queryId)) {
            this.pendingQueries.get(message.queryId)!(message.state);
            this.pendingQueries.delete(message.queryId);
        }
    }

    query(): Promise<SystemState> {
        return new Promise(resolve => {
            const queryId = this.queryIdCounter++;
            this.pendingQueries.set(queryId, resolve);
            const msg: QueryMessage = { type: 'query', queryId };
            void this.connection.postMessage(msg);
        });
    }

    tick(deltaTime: number, input: InputState): Promise<void> {
        const msg: ClientToServerMessage = { type: 'tick', deltaTime, input };
        void this.connection.postMessage(msg);
        return Promise.resolve();
    }

    setTimeScale(scale: number): Promise<void> {
        const msg: ClientToServerMessage = { type: 'setTimeScale', scale };
        return this.connection.postMessage(msg);
    }

    addBody(body: Omit<Body, 'id'>): Promise<void> {
        const msg: ClientToServerMessage = { type: 'addBody', body };
        return this.connection.postMessage(msg);
    }

    autoLand(targetBodyId: string | null): Promise<void> {
        const msg: ClientToServerMessage = { type: 'autoLand', targetBodyId };
        return this.connection.postMessage(msg);
    }

    teleportToSurface(targetBodyId: string | null): Promise<void> {
        const msg: ClientToServerMessage = { type: 'teleportToSurface', targetBodyId };
        return this.connection.postMessage(msg);
    }
}



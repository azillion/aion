import type { Authority, InputState } from '@shared/authority';
import type { Body, SystemState } from '@shared/types';

// This interface defines the contract for any transport layer (in-memory, websocket, etc.)
export interface IAuthorityConnection {
    postMessage(message: any): Promise<void>;
    onMessage: ((message: any) => void) | null;
}

export class ClientAuthority implements Authority {
    private connection: IAuthorityConnection;
    private pendingQueries: Map<number, (value: SystemState) => void> = new Map();
    private queryIdCounter = 0;

    constructor(connection: IAuthorityConnection) {
        this.connection = connection;
        this.connection.onMessage = this.handleMessage.bind(this);
    }

    private handleMessage(message: any): void {
        if (message.type === 'queryResult' && this.pendingQueries.has(message.queryId)) {
            this.pendingQueries.get(message.queryId)!(message.state);
            this.pendingQueries.delete(message.queryId);
        }
    }

    query(): Promise<SystemState> {
        return new Promise(resolve => {
            const queryId = this.queryIdCounter++;
            this.pendingQueries.set(queryId, resolve);
            void this.connection.postMessage({ type: 'query', queryId });
        });
    }

    tick(deltaTime: number, input: InputState): Promise<void> {
        void this.connection.postMessage({ type: 'tick', deltaTime, input });
        return Promise.resolve();
    }

    setTimeScale(scale: number): void {
        void this.connection.postMessage({ type: 'setTimeScale', scale });
    }

    addBody(body: Omit<Body, 'id'>): void {
        void this.connection.postMessage({ type: 'addBody', body });
    }

    autoLand(targetBodyId: string | null): void {
        void this.connection.postMessage({ type: 'autoLand', targetBodyId });
    }

    teleportToSurface(targetBodyId: string | null): void {
        void this.connection.postMessage({ type: 'teleportToSurface', targetBodyId });
    }
}



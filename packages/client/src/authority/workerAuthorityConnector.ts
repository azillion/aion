import type { IAuthorityConnection } from './clientAuthority';
import type { ClientToServerMessage, ServerToClientMessage } from '@shared/messages';

export class WorkerAuthorityConnector implements IAuthorityConnection {
    private worker: Worker;
    public onMessage: ((message: ServerToClientMessage) => void) | null = null;

    constructor() {
        // The { type: 'module' } is crucial for allowing the worker to use import statements.
        this.worker = new Worker(new URL('../../../server/src/worker.ts', import.meta.url), {
            type: 'module',
        });

        this.worker.onmessage = (e: MessageEvent<ServerToClientMessage>) => {
            this.onMessage?.(e.data);
        };
    }

    async postMessage(message: ClientToServerMessage): Promise<void> {
        this.worker.postMessage(message);
    }
}



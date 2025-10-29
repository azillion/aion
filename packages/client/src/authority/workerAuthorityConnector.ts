import type { IAuthorityConnection } from './clientAuthority';

export class WorkerAuthorityConnector implements IAuthorityConnection {
    private worker: Worker;
    public onMessage: ((message: any) => void) | null = null;

    constructor() {
        // The { type: 'module' } is crucial for allowing the worker to use import statements.
        this.worker = new Worker(new URL('../../../server/src/worker.ts', import.meta.url), {
            type: 'module',
        });

        this.worker.onmessage = (e: MessageEvent) => {
            this.onMessage?.(e.data);
        };
    }

    async postMessage(message: any): Promise<void> {
        this.worker.postMessage(message);
    }
}



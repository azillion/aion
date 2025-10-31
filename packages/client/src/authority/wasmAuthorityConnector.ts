import type { IAuthorityConnection } from './provider';
import type { ClientToServerMessage, ServerToClientMessage } from '@shared/messages';

export class WasmAuthorityConnector implements IAuthorityConnection {
    private worker: Worker;
    public onMessage: ((message: ServerToClientMessage) => void) | null = null;
    private isReady = false;
    private messageQueue: ClientToServerMessage[] = [];

    constructor() {
        this.worker = new Worker(new URL('./wasm.worker.ts', import.meta.url), { type: 'module' });
        this.worker.onmessage = (e: MessageEvent<any>) => {
            const data = e.data as any;
            if (data && data.type === 'workerReady') {
                this.isReady = true;
                // flush queued messages
                for (const msg of this.messageQueue) {
                    this.worker.postMessage(msg);
                }
                this.messageQueue = [];
                return;
            }
            this.onMessage?.(data as ServerToClientMessage);
        };
    }

    async postMessage(message: ClientToServerMessage): Promise<void> {
        if (this.isReady) {
            this.worker.postMessage(message);
        } else {
            this.messageQueue.push(message);
        }
    }
}



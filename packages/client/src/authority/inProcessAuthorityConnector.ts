import type { LocalAuthority } from '@server/local';
import type { IAuthorityConnection } from './clientAuthority';
import type { ClientToServerMessage, ServerToClientMessage } from '@shared/messages';

export class InProcessAuthorityConnector implements IAuthorityConnection {
    private server: LocalAuthority;
    public onMessage: ((message: ServerToClientMessage) => void) | null = null;

    constructor(server: LocalAuthority) {
        this.server = server;
    }

    async postMessage(message: ClientToServerMessage): Promise<void> {
        // Simulate the async nature of a network call
        await Promise.resolve();

        switch (message.type) {
            case 'query': {
                const state = await this.server.query();
                this.onMessage?.({ type: 'queryResult', queryId: message.queryId, state });
                break;
            }
            case 'tick':
                await this.server.tick(message.deltaTime, message.input);
                break;
            case 'setTimeScale':
                await this.server.setTimeScale(message.scale);
                break;
            case 'addBody':
                await this.server.addBody(message.body);
                break;
            case 'autoLand':
                await this.server.autoLand(message.targetBodyId);
                break;
            case 'teleportToSurface':
                await this.server.teleportToSurface(message.targetBodyId);
                break;
        }
    }
}



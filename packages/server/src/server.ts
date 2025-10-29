import { LocalAuthority } from './local';
import type { ClientToServerMessage, QueryResultMessage } from '@shared/messages';

export function runServer(postResponse: (message: any) => void) {
    const server = new LocalAuthority();

    const handleMessage = async (message: ClientToServerMessage) => {
        switch (message.type) {
            case 'query': {
                const state = await server.query();
                const response: QueryResultMessage = { type: 'queryResult', queryId: message.queryId, state };
                postResponse(response);
                break;
            }
            case 'tick':
                await server.tick(message.deltaTime, message.input);
                break;
            case 'setTimeScale':
                await server.setTimeScale(message.scale);
                break;
            case 'addBody':
                await server.addBody(message.body);
                break;
            case 'autoLand':
                await server.autoLand(message.targetBodyId);
                break;
            case 'teleportToSurface':
                await server.teleportToSurface(message.targetBodyId);
                break;
        }
    };

    return handleMessage;
}



import { LocalAuthority } from './local';
import type { ClientToServerMessage, QueryResultMessage } from '@shared/messages';

const server = new LocalAuthority();

self.onmessage = async (e: MessageEvent<ClientToServerMessage>) => {
  const message = e.data;
  switch (message.type) {
    case 'query': {
      const state = await server.query();
      const response: QueryResultMessage = { type: 'queryResult', queryId: message.queryId, state };
      self.postMessage(response);
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
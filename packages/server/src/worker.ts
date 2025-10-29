import { LocalAuthority } from './local';

const server = new LocalAuthority();

self.onmessage = async (e: MessageEvent) => {
  const message = e.data;
  switch (message.type) {
    case 'query': {
      const state = await server.query();
      self.postMessage({ type: 'queryResult', queryId: message.queryId, state });
      break;
    }
    case 'tick':
      await server.tick(message.deltaTime, message.input);
      break;
    case 'setTimeScale':
      server.setTimeScale(message.scale);
      break;
    case 'addBody':
      server.addBody(message.body);
      break;
    case 'autoLand':
      server.autoLand(message.targetBodyId);
      break;
    case 'teleportToSurface':
      server.teleportToSurface(message.targetBodyId);
      break;
  }
};



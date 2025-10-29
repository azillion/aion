import { runServer } from './server';

const handleMessage = runServer(self.postMessage.bind(self));

self.onmessage = (e: MessageEvent) => {
  handleMessage(e.data);
};
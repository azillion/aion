import { WebSocketServer, WebSocket } from 'ws';
import { FfiSimulationService } from './ffiSimulationService';
import type { ISimulationService } from './simulationService';

async function main() {
    const service: ISimulationService = new FfiSimulationService();
    await service.initialize();

    const wss = new WebSocketServer({ port: 9001 });
    // eslint-disable-next-line no-console
    console.log('WebSocket server listening on ws://localhost:9001');

    const tick_rate_ms = 1000 / 10;
    const mainLoop = setInterval(() => {
        service.tick();
        const data = service.getLatestData();
        if (data) {
            wss.clients.forEach(client => {
                if (client.readyState === WebSocket.OPEN) {
                    client.send(data);
                }
            });
        }
    }, tick_rate_ms);

    process.on('SIGINT', () => {
        // eslint-disable-next-line no-console
        console.log('Caught interrupt signal, shutting down.');
        clearInterval(mainLoop);
        service.shutdown();
        wss.close();
        process.exit(0);
    });
}

main().catch(err => {
    // eslint-disable-next-line no-console
    console.error('Server failed to start:', err);
    process.exit(1);
});



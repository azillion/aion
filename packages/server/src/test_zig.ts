import { FfiSimulationService } from './ffiSimulationService';

async function main() {
    // eslint-disable-next-line no-console
    console.log('Starting Zig pointer debug harness...');
    const svc = new FfiSimulationService();
    await svc.initialize();

    // One tick to advance the simulation
    svc.tick();

    // This call should trigger Zig's debug prints and then exit(0)
    // eslint-disable-next-line no-console
    console.log('About to call getLatestData (Zig will exit after logging)...');
    // Intentionally ignore the returned data; process should exit from Zig
    void svc.getLatestData();
}

main().catch(err => {
    // eslint-disable-next-line no-console
    console.error('Debug harness failed:', err);
    process.exit(1);
});



import type { IAuthorityConnection } from './clientAuthority';
import { WorkerAuthorityConnector } from './workerAuthorityConnector';

export interface IAuthorityProvider {
    createConnection(): Promise<IAuthorityConnection>;
}

export class WebWorkerAuthorityProvider implements IAuthorityProvider {
    public async createConnection(): Promise<IAuthorityConnection> {
        return Promise.resolve(new WorkerAuthorityConnector());
    }
}



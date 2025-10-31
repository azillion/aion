import type { IAuthorityConnection } from './clientAuthority';
import { WasmAuthorityConnector } from './wasmAuthorityConnector';

export interface IAuthorityProvider {
    createConnection(): Promise<IAuthorityConnection>;
}

export class WasmAuthorityProvider implements IAuthorityProvider {
    async createConnection(): Promise<IAuthorityConnection> {
        return Promise.resolve(new WasmAuthorityConnector());
    }
}



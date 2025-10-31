import type { IAuthorityConnection, IAuthorityProvider } from './provider';
import { WasmAuthorityConnector } from './wasmAuthorityConnector';

export class WasmAuthorityProvider implements IAuthorityProvider {
    async createConnection(): Promise<IAuthorityConnection> {
        return Promise.resolve(new WasmAuthorityConnector());
    }
}



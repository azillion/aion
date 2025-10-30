/**
 * Defines the contract for a simulation service.
 * The main server loop will interact with this interface,
 * remaining agnostic to the underlying implementation (FFI, Redis, etc.).
 */
export interface ISimulationService {
  /**
   * Performs one-time asynchronous initialization.
   * @returns A promise that resolves when the service is ready.
   */
  initialize(): Promise<void>;

  /**
   * Advances the simulation by one discrete step.
   * This should be a fast, non-blocking operation.
   */
  tick(): void;

  /**
   * Retrieves the latest complete simulation state as a binary buffer.
   * This may be a blocking operation that waits for the GPU.
   * @returns A Buffer containing the raw simulation data, or null if not ready.
   */
  getLatestData(): Buffer | null;

  /**
   * Performs graceful shutdown and releases all resources.
   */
  shutdown(): void;
}



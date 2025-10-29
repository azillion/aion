export interface ISimulationHost {
  initialize(device: GPUDevice, queue: GPUQueue): Promise<void>;
  tick(): void;
  getOutputTextureView(): GPUTextureView;
}

export type SimulationHandle = number;



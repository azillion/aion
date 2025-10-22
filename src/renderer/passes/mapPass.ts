import type { WebGPUCore } from '../core';
import type { Scene } from '../scene';
import type { IRenderPass, RenderContext } from '../types';
import type { Theme } from '../../shared/types';
import mapShaderWGSL from '../shaders/map.wgsl?raw';
import cameraWGSL from '../shaders/camera.wgsl?raw';
import backgroundShaderWGSL from '../shaders/background.wgsl?raw';

export class MapPass implements IRenderPass {
  private pipeline!: GPURenderPipeline;
  private backgroundPipeline!: GPURenderPipeline;

  public initialize(core: WebGPUCore, _scene: Scene): void {
    const module = core.device.createShaderModule({ code: cameraWGSL + mapShaderWGSL });
    this.pipeline = core.device.createRenderPipeline({
      label: 'Map Pipeline',
      layout: 'auto',
      vertex: { module, entryPoint: 'vertexMain' },
      fragment: { 
        module, 
        entryPoint: 'fragmentMain', 
        targets: [{ 
          format: core.presentationFormat,
          blend: {
            color: { srcFactor: 'src-alpha', dstFactor: 'one-minus-src-alpha', operation: 'add' },
            alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add' },
          }
        }] 
      },
      primitive: { topology: 'triangle-strip' },
    });

    // Background pipeline to draw the precomputed orbits texture
    const bgModule = core.device.createShaderModule({ code: backgroundShaderWGSL });
    this.backgroundPipeline = core.device.createRenderPipeline({
        label: 'Map Background Pipeline',
        layout: 'auto',
        vertex: { module: bgModule, entryPoint: 'vertexMain' },
        fragment: {
            module: bgModule,
            entryPoint: 'fragmentMain',
            targets: [{ format: core.presentationFormat }],
        },
    });
  }

  public run(encoder: GPUCommandEncoder, context: RenderContext, _theme: Theme, showOrbits: boolean): void {
    const { core, scene, orbitsTexture } = context;

    // 1. Shared camera buffer is already populated by Scene.update

    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: core.context.getCurrentTexture().createView(),
        loadOp: 'clear',
        storeOp: 'store',
        clearValue: { r: 0.0, g: 0.0, b: 0.0, a: 1.0 }
      }],
    });
    
    // 2. Draw the Orbits Texture as the background
    if (showOrbits) {
        const sampler = core.device.createSampler();
        const backgroundBindGroup = core.device.createBindGroup({
            label: 'Map Background BindGroup',
            layout: this.backgroundPipeline.getBindGroupLayout(0),
            entries: [
                { binding: 0, resource: sampler },
                { binding: 1, resource: orbitsTexture.createView() }
            ]
        });
        pass.setPipeline(this.backgroundPipeline);
        pass.setBindGroup(0, backgroundBindGroup);
        pass.draw(6);
    }

    // 3. Draw the planet icons on top
    const planetBindGroup = core.device.createBindGroup({
      label: 'Map Planet Bind Group',
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: scene.sharedCameraUniformBuffer } },
        { binding: 1, resource: { buffer: scene.spheresBuffer } },
      ],
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, planetBindGroup);
    pass.draw(4, scene.lastKnownBodyCount);
    pass.end();
  }
}



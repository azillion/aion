import { defineConfig } from 'vite';
import { viteStaticCopy } from 'vite-plugin-static-copy';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const rootDir = path.dirname(fileURLToPath(import.meta.url));

export default defineConfig({
  worker: {
    format: 'es',
  },
  resolve: {
    alias: {
      '@shared': path.resolve(rootDir, 'packages/shared/src'),
      '@client': path.resolve(rootDir, 'packages/client/src'),
      '@server': path.resolve(rootDir, 'packages/server/src'),
    },
  },
  server: {
    fs: {
      allow: ['..'],
    },
  },
  plugins: [
    viteStaticCopy({
      targets: [
        {
          src: 'packages/game-sim/zig-out/bin/game-sim.wasm',
          dest: 'packages/game-sim/zig-out/bin',
        },
      ],
    }),
  ],
});



import { defineConfig } from 'vite';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const rootDir = path.dirname(fileURLToPath(import.meta.url));

export default defineConfig({
  resolve: {
    alias: {
      '@shared': path.resolve(rootDir, 'packages/shared/src'),
      '@client': path.resolve(rootDir, 'packages/client/src'),
      '@server': path.resolve(rootDir, 'packages/server/src'),
    },
  },
});



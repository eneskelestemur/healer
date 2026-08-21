import { defineConfig, loadEnv } from 'vite'
import react from '@vitejs/plugin-react'
import path from 'path'

// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, path.resolve(process.cwd(), '..'), '');
  return {
    plugins: [react()],
    envDir: '../',
    build: {
      // Output compiled assets directly to the python package directory.
      // The directory also holds the checked-in ketcher bundle, so it is
      // neither emptied nor overwritten from public/ (which only symlinks it).
      outDir: '../healer/web/static',
      emptyOutDir: false,
      copyPublicDir: false,
    },
    server: {
      host: true, // Needed for Docker
      proxy: {
        // Proxy API requests to the Python Backend during development
        '/api': {
          target: env.VITE_API_TARGET || 'http://127.0.0.1:8000',
          changeOrigin: true,
          secure: false,
        }
      }
    }
  }
})

/// <reference types="vite/client" />

interface ImportMetaEnv {
  readonly VITE_MAX_DISPLAY_ROWS?: string;
  // add other env variables here...
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}

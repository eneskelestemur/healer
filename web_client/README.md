# HEALER Web Client

This directory contains the React source code for the HEALER web interface.

## Ketcher

`public/ketcher` is a symlink to `healer/web/static/ketcher`, the copy shipped in
the Python package. Keeping one copy means the editor is served from the same
files in development and after a pip install.

On Windows, symlinks need Developer Mode or `git config core.symlinks true`. If
the link does not resolve, copy `healer/web/static/ketcher` to
`web_client/public/ketcher` instead.

## Running Locally (without Docker)

### Local Mode (Simple)

Only 2 terminals needed - no Redis/Celery required:

```bash
# Terminal 1: Backend
cd <repo root>
python -m uvicorn healer.web.app:app --reload --port 8000

# Terminal 2: Frontend
cd web_client
npm run dev
```

### Server Mode (with Celery/Redis)

4 terminals needed for full async job support:

| Terminal | Command | Purpose |
|----------|---------|---------|
| 1 | `redis-server` | Message broker & result backend |
| 2 | `celery -A healer.web.celery_worker worker --loglevel=info` | Async job worker |
| 3 | `HEALER_SERVER_MODE=true uvicorn healer.web.app:app --reload` | REST API |
| 4 | `npm run dev` | React frontend |

```bash
# Terminal 1: Redis
redis-server

# Terminal 2: Celery Worker
cd <repo root>
celery -A healer.web.celery_worker worker --loglevel=info

# Terminal 3: Backend (with Celery enabled)
cd <repo root>
HEALER_SERVER_MODE=true python -m uvicorn healer.web.app:app --reload --port 8000

# Terminal 4: Frontend
cd web_client
npm run dev
```

Then open **http://localhost:5173** in your browser.

> **Note**: The `HEALER_SERVER_MODE=true` environment variable is required for server mode. Without it, the backend runs in local (synchronous) mode.

## Development

1. Install dependencies:
   ```bash
   npm install
   ```

2. Start the development server:
   ```bash
   npm run dev
   ```

## Building for Production

To build the frontend for inclusion in the Python package:

```bash
npm run build
```

This should be configured to output files to `../healer/web/static`.

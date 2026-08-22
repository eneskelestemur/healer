# Web Interface

An optional browser interface with a structure editor, for interactive use.

```bash
pip install 'mol-healer[web]'
healer-ui
```

Then open <http://localhost:8000>.

| Option | Description | Default |
|--------|-------------|---------|
| `--host` | Bind address | `0.0.0.0` |
| `--port` | Port | `8000` |

The frontend is bundled with the package, so a pip install is all that is needed
— there is nothing to build.

## Local mode vs. server mode

By default jobs run in the same process and are held in memory, which suits a
single user on a laptop.

Setting `HEALER_SERVER_MODE=true` switches to Celery with Redis, so jobs are
queued and executed by workers. Server mode also enforces the caps below, so a
shared deployment cannot be overwhelmed by one request.

| Variable | Caps | Default |
|----------|------|---------|
| `HEALER_LIMIT_MAX_EVALS` | Reaction attempts per composition | 1000000 |
| `HEALER_LIMIT_MAX_PRODUCTS` | Products per composition | 1000 |
| `HEALER_LIMIT_MAX_TOTAL` | Products in total | 50000 |
| `HEALER_LIMIT_SIM_MIN` / `_MAX` | Similarity threshold range | 0.65 / 1.0 |
| `HEALER_LIMIT_MAX_BBS` | Building blocks per fragment | 10 |
| `HEALER_LIMIT_N_COMP` | Compositions | 50 |
| `HEALER_LIMIT_RETRO_DEPTH` | Retrosynthesis depth | 2 |
| `HEALER_LIMIT_MIN_FRAG` | Minimum fragment size | 7 |
| `HEALER_LIMIT_MAX_RXN_TAGS` | Reaction tags per request | 10 |

Requests exceeding a cap are clamped rather than rejected.

## Docker

`docker-compose.yml` brings up the API, a Redis broker, and a Celery worker in
server mode:

```bash
docker compose up
```

Copy `.env.example` to `.env` first to set the building block directory and
limits.

## HTTP API

Enumeration is asynchronous: submit a job, poll it, then download the results.
Every path below is served under the `/api` prefix.

| Method | Path | Purpose |
|--------|------|---------|
| `POST` | `/enumerate/molecule` | Submit a molecule job |
| `POST` | `/enumerate/site` | Submit a site job |
| `GET` | `/jobs/{job_id}` | Poll status and results |
| `POST` | `/jobs/{job_id}/cancel` | Cancel a job |
| `GET` | `/jobs/{job_id}/download` | Download results as CSV |
| `GET` | `/health` | Liveness check |
| `GET` | `/info/mode` | `local` or `celery` |
| `GET` | `/info/limits` | Active server limits |
| `GET` | `/info/building-blocks` | Available building block sources |
| `GET` | `/utils/reaction-tags` | Available reaction tags |
| `POST` | `/utils/smiles-to-mol` | SMILES to molfile, for the editor |
| `POST` | `/utils/render-mol-with-indices` | SVG with atom indices |
| `POST` | `/utils/render-result` | SVG of a product with its building blocks |

Interactive API documentation is served at `/docs`.

```bash
JOB=$(curl -s -X POST localhost:8000/api/enumerate/molecule \
  -H 'Content-Type: application/json' \
  -d '{"molecule": "CC(=O)Nc1ccccc1", "bb_source": "test"}' | jq -r .job_id)

curl -s localhost:8000/api/jobs/$JOB | jq .status
```

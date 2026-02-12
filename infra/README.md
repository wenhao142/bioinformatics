## Run local stack (T0.1 scaffold)

```bash
cd infra
cp .env.example .env  # optional; compose has sensible defaults
docker compose up --build
```

Services (host ports remapped to avoid conflicts):
- web: http://localhost:13000 (Next.js)
- api: http://localhost:18000/health (FastAPI)
- postgres: localhost:15432 (postgres/postgres, db=omics)
- minio: http://localhost:19001 (minioadmin/minioadmin123)

Expected: UI loads, API health returns `{"status":"ok"}`.

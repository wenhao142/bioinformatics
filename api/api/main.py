from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from api import auth, rbac, audit, datasets, variants

app = FastAPI(title="AD Locus Evidence API")

# Allow local dev cross-origin (web runs on 13000, API on 18000)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(auth.router)
app.include_router(rbac.router)
app.include_router(audit.router)
app.include_router(datasets.router)
app.include_router(variants.router)


@app.get("/health")
def health():
    return {"status": "ok"}

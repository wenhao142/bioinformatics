import os
import gzip
import psycopg2
from io import BytesIO
from typing import List, TypedDict, Optional

import jwt
from fastapi import APIRouter, Depends, File, HTTPException, Query, UploadFile, Response, Security
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

from api.auth import current_user, JWT_SECRET

router = APIRouter(prefix="/variants", tags=["variants"])
bearer = HTTPBearer(auto_error=False)


class Variant(TypedDict):
    chr: str
    pos: int
    ref: str
    alt: str
    qual: Optional[float]
    filter: str
    gene: Optional[str]


class PgVariantStore:
    def __init__(self):
        self.conn = psycopg2.connect(self._dsn())
        self.conn.autocommit = True
        self._ensure_table()

    def _dsn(self) -> str:
        dsn = os.getenv("DATABASE_URL")
        if dsn:
            return dsn
        return (
            f"dbname={os.getenv('POSTGRES_DB', 'omics')} "
            f"user={os.getenv('POSTGRES_USER', 'postgres')} "
            f"password={os.getenv('POSTGRES_PASSWORD', 'postgres')} "
            f"host={os.getenv('POSTGRES_HOST', 'localhost')} "
            f"port={os.getenv('POSTGRES_PORT', '5432')}"
        )

    def _ensure_table(self):
        with self.conn.cursor() as cur:
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS variants (
                    id SERIAL PRIMARY KEY,
                    chr TEXT NOT NULL,
                    pos INTEGER NOT NULL,
                    ref TEXT NOT NULL,
                    alt TEXT NOT NULL,
                    qual DOUBLE PRECISION,
                    filter_status TEXT
                );
                CREATE INDEX IF NOT EXISTS variants_chr_pos_idx ON variants(chr, pos);
                """
            )

    def insert_many(self, rows: List[Variant]):
        with self.conn.cursor() as cur:
            cur.executemany(
                "INSERT INTO variants (chr, pos, ref, alt, qual, filter_status) VALUES (%s,%s,%s,%s,%s,%s)",
                [
                    (v["chr"], v["pos"], v["ref"], v["alt"], v["qual"], v["filter"])
                    for v in rows
                ],
            )

    def query(self, chr: str, start: int, end: int) -> List[Variant]:
        with self.conn.cursor() as cur:
            cur.execute(
                "SELECT chr, pos, ref, alt, qual, filter_status FROM variants WHERE chr=%s AND pos BETWEEN %s AND %s ORDER BY pos LIMIT 1000",
                (chr, start, end),
            )
            return [
                {
                    "chr": r[0],
                    "pos": r[1],
                    "ref": r[2],
                    "alt": r[3],
                    "qual": r[4],
                    "filter": r[5],
                }
                for r in cur.fetchall()
            ]

    def stats(self):
        with self.conn.cursor() as cur:
            cur.execute("SELECT count(*) FROM variants")
            total = cur.fetchone()[0]
        return {"total_variants": total}

    def all_variants(self) -> List[Variant]:
        with self.conn.cursor() as cur:
            cur.execute("SELECT chr, pos, ref, alt, qual, filter_status FROM variants LIMIT 100000")
            return [
                {
                    "chr": r[0],
                    "pos": r[1],
                    "ref": r[2],
                    "alt": r[3],
                    "qual": r[4],
                    "filter": r[5],
                    "gene": None,
                }
                for r in cur.fetchall()
            ]


class MemoryVariantStore:
    def __init__(self):
        self.rows: List[Variant] = []

    def insert_many(self, rows: List[Variant]):
        self.rows.extend(rows)

    def query(self, chr: str, start: int, end: int) -> List[Variant]:
        return [v for v in self.rows if v["chr"] == chr and start <= v["pos"] <= end][:1000]

    def stats(self):
        return {"total_variants": len(self.rows)}

    def all_variants(self) -> List[Variant]:
        return list(self.rows)


def get_store():
    if os.getenv("VARIANTS_USE_DB", "false").lower() in ("1", "true", "yes", "on"):
        try:
            return PgVariantStore()
        except Exception as exc:  # pragma: no cover - fallback
            raise RuntimeError(f"Failed to init Postgres variant store: {exc}")
    return MemoryVariantStore()


STORE = get_store()


def parse_vcf_bytes(data: bytes) -> List[Variant]:
    rows: List[Variant] = []
    buf = BytesIO(data)
    opener = gzip.open if data[:2] == b"\x1f\x8b" else lambda x: x
    fh = opener(buf)
    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 6:
            continue
        chrom, pos, _id, ref, alt, qual, filt = parts[:7]
        try:
            pos_int = int(pos)
        except ValueError:
            continue
        try:
            qual_val = float(qual) if qual not in (".", "") else None
        except ValueError:
            qual_val = None
        for a in alt.split(","):
            rows.append(
                {
                    "chr": chrom,
                    "pos": pos_int,
                    "ref": ref,
                    "alt": a,
                    "qual": qual_val,
                    "filter": filt,
                    "gene": None,
                }
            )
    return rows


# Minimal built-in gene intervals for MVP nearest-gene mapping
GeneInterval = TypedDict("GeneInterval", {"chr": str, "start": int, "end": int, "gene": str})
GENE_INTERVALS: List[GeneInterval] = [
    {"chr": "chr1", "start": 1, "end": 500, "gene": "GENE1"},
    {"chr": "chr1", "start": 600, "end": 900, "gene": "GENE2"},
    {"chr": "chr2", "start": 100, "end": 300, "gene": "GENE3"},
]


def nearest_gene(chr: str, pos: int) -> Optional[str]:
    candidates = [g for g in GENE_INTERVALS if g["chr"] == chr]
    if not candidates:
        return None
    best_gene = None
    best_dist = 10**9
    for g in candidates:
        if g["start"] <= pos <= g["end"]:
            return g["gene"]
        dist = min(abs(pos - g["start"]), abs(pos - g["end"]))
        if dist < best_dist:
            best_dist = dist
            best_gene = g["gene"]
    return best_gene


@router.post("/ingest")
def ingest_vcf(file: UploadFile = File(...), project_id: str | None = None, user=Depends(current_user)):
    data = file.file.read()
    variants = parse_vcf_bytes(data)
    if not variants:
        raise HTTPException(status_code=400, detail="No variants parsed")
    STORE.insert_many(variants)
    return {"ingested": len(variants), "project_id": project_id}


@router.get("")
def list_variants(
    chr: str = Query(..., description="chromosome, e.g., chr1"),
    start: int = Query(1),
    end: int = Query(250_000_000),
    annotate: bool = Query(False, description="attach nearest gene"),
    user=Depends(current_user),
):
    variants = STORE.query(chr, start, end)
    if annotate:
        for v in variants:
            v["gene"] = nearest_gene(v["chr"], v["pos"])
    return {"variants": variants}


@router.get("/stats")
def variant_stats(user=Depends(current_user)):
    return STORE.stats()


@router.get("/bed", response_class=Response)
def variants_bed(
    chr: str = Query(..., description="chromosome, e.g., chr1"),
    start: int = Query(1),
    end: int = Query(250_000_000),
    token: Optional[str] = Query(None, description="JWT token (optional if Authorization header present)"),
    credentials: HTTPAuthorizationCredentials | None = Security(bearer),
):
    # authorize via header or token query
    raw_token = token
    if credentials and credentials.credentials:
        raw_token = credentials.credentials
    if not raw_token:
        raise HTTPException(status_code=401, detail="Missing token")
    try:
        jwt.decode(raw_token, JWT_SECRET, algorithms=["HS256"])
    except jwt.PyJWTError:
        raise HTTPException(status_code=401, detail="Invalid token")

    variants = STORE.query(chr, start, end)
    lines = []
    for v in variants:
        bed_start = max(0, v["pos"] - 1)
        bed_end = v["pos"]
        name = f"{v['ref']}>{v['alt']}"
        lines.append(f"{v['chr']}\t{bed_start}\t{bed_end}\t{name}\t0\t.")
    body = "\n".join(lines) + ("\n" if lines else "")
    return Response(content=body, media_type="text/plain")


@router.get("/nearest")
def nearest(
    chr: str = Query(..., description="chromosome, e.g., chr1"),
    pos: int = Query(...),
    user=Depends(current_user),
):
    return {"gene": nearest_gene(chr, pos)}

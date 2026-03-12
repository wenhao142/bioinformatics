import hashlib
import gzip
import os
import tempfile
import time
from pathlib import Path
from typing import List, TypedDict

from fastapi import APIRouter, Depends, File, HTTPException, UploadFile

from api.auth import current_user

try:
    from minio import Minio
    from minio.error import S3Error
except Exception:  # pragma: no cover - minio optional in tests
    Minio = None
    S3Error = Exception


class DatasetRecord(TypedDict):
    id: int
    filename: str
    sha256: str
    size: int
    uri: str
    uploaded_by: str
    uploaded_at: float
    project_id: str | None
    input_role: str
    canonical_type: str | None
    validation_status: str
    validation_detail: str


router = APIRouter(prefix="/datasets", tags=["datasets"])

_DATASETS: List[DatasetRecord] = []
_NEXT_ID = 1


def _compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_text_preview(path: Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return [handle.readline().rstrip("\n") for _ in range(4)]
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        return [handle.readline().rstrip("\n") for _ in range(4)]


def _infer_dataset_type(path: Path) -> tuple[str, str | None]:
    name = path.name.lower()
    if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "reads", "reads.fastq.gz"
    if name.endswith((".fasta", ".fa", ".fna")):
        return "reference", "reference.genome.fasta"
    if name.endswith((".gtf", ".gff", ".gff3")):
        return "annotation", "annotation.gtf"
    if name.endswith(".bam"):
        return "alignment", "align.bam"
    if name.endswith(".vcf.gz"):
        return "variants", "variants.vcf.gz"
    if name.endswith(".tsv"):
        return "table", "expression.diff_table.tsv"
    return "unknown", None


def _validate_uploaded_content(path: Path, input_role: str) -> tuple[str, str]:
    try:
        lines = _read_text_preview(path)
    except Exception as exc:
        if input_role == "alignment":
            return "validated", "Binary alignment accepted by extension."
        return "invalid", f"Cannot inspect file content: {exc}"

    if input_role == "reads":
        if len(lines) >= 4 and lines[0].startswith("@") and lines[2].startswith("+"):
            return "validated", "FASTQ structure looks valid."
        return "invalid", "Expected FASTQ structure (@ header and + quality separator)."
    if input_role == "reference":
        if lines and lines[0].startswith(">"):
            return "validated", "FASTA structure looks valid."
        return "invalid", "Expected FASTA content starting with '>'."
    if input_role == "annotation":
        content_lines = [line for line in lines if line and not line.startswith("#")]
        if content_lines and len(content_lines[0].split("\t")) >= 8:
            return "validated", "Annotation table looks valid."
        return "invalid", "Expected GTF/GFF tab-delimited annotation content."
    if input_role == "table":
        if any("\t" in line for line in lines if line):
            return "validated", "Tabular text accepted."
        return "invalid", "Expected tab-delimited table content."
    if input_role == "variants":
        if lines and (lines[0].startswith("##") or lines[0].startswith("#CHROM")):
            return "validated", "VCF header detected."
        return "invalid", "Expected VCF header content."
    if input_role == "alignment":
        return "validated", "Alignment accepted by extension."
    return "invalid", "Unsupported raw file type for workflow input."


def _get_minio_client():
    if Minio is None:
        return None
    endpoint = os.getenv("MINIO_ENDPOINT", "http://localhost:9000")
    access_key = os.getenv("MINIO_ROOT_USER", "minioadmin")
    secret_key = os.getenv("MINIO_ROOT_PASSWORD", "minioadmin123")
    secure = endpoint.startswith("https")
    # strip scheme for Minio client
    endpoint_clean = endpoint.replace("https://", "").replace("http://", "")
    return Minio(endpoint_clean, access_key=access_key, secret_key=secret_key, secure=secure)


def _ensure_bucket(client: Minio, bucket: str):
    if client.bucket_exists(bucket):
        return
    client.make_bucket(bucket)


def _save_to_minio(tmp_path: Path, object_name: str) -> str:
    client = _get_minio_client()
    if client is None:
        raise RuntimeError("MinIO client not available")
    bucket = os.getenv("DATASET_BUCKET", "datasets")
    _ensure_bucket(client, bucket)
    client.fput_object(bucket, object_name, str(tmp_path))
    return f"s3://{bucket}/{object_name}"


def _save_to_local(tmp_path: Path, object_name: str) -> str:
    base = Path(os.getenv("LOCAL_DATASET_DIR", "/tmp/datasets"))
    base.mkdir(parents=True, exist_ok=True)
    dest = base / object_name
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(tmp_path.read_bytes())
    return str(dest)


def _save_file(upload: UploadFile, user_email: str, project_id: str | None) -> DatasetRecord:
    global _NEXT_ID
    use_minio = os.getenv("USE_MINIO", "true").lower() in ("1", "true", "yes", "on")
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp.write(upload.file.read())
        tmp_path = Path(tmp.name)
    sha = _compute_sha256(tmp_path)
    size = tmp_path.stat().st_size
    input_role, canonical_type = _infer_dataset_type(Path(upload.filename or ""))
    validation_status, validation_detail = _validate_uploaded_content(tmp_path, input_role)
    if input_role == "unknown" or canonical_type is None:
        tmp_path.unlink(missing_ok=True)
        raise HTTPException(status_code=400, detail="Unsupported raw file type for workflow input")
    if validation_status != "validated":
        tmp_path.unlink(missing_ok=True)
        raise HTTPException(status_code=400, detail=validation_detail)
    object_name = f"{sha}/{upload.filename}"
    if use_minio:
        try:
            uri = _save_to_minio(tmp_path, object_name)
        except Exception as exc:
            tmp_path.unlink(missing_ok=True)
            raise HTTPException(status_code=500, detail=f"MinIO upload failed: {exc}")
    else:
        uri = _save_to_local(tmp_path, object_name)
    tmp_path.unlink(missing_ok=True)

    record: DatasetRecord = {
        "id": _NEXT_ID,
        "filename": upload.filename,
        "sha256": sha,
        "size": size,
        "uri": uri,
        "uploaded_by": user_email,
        "uploaded_at": time.time(),
        "project_id": project_id,
        "input_role": input_role,
        "canonical_type": canonical_type,
        "validation_status": validation_status,
        "validation_detail": validation_detail,
    }
    _DATASETS.append(record)
    _NEXT_ID += 1
    return record


@router.post("/upload")
def upload_dataset(
    file: UploadFile = File(...),
    project_id: str | None = None,
    user=Depends(current_user),
):
    rec = _save_file(file, user["email"], project_id)
    return {"dataset": rec}


@router.get("")
def list_datasets(user=Depends(current_user)):
    if user["role"] == "admin":
        return {"datasets": _DATASETS}
    return {"datasets": [d for d in _DATASETS if d["uploaded_by"] == user["email"]]}

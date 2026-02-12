import hashlib
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


router = APIRouter(prefix="/datasets", tags=["datasets"])

_DATASETS: List[DatasetRecord] = []
_NEXT_ID = 1


def _compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


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

"""
Poll a MinIO bucket for new VCF/VCF.GZ files, then:
- downloads to /tmp
- extracts core fields via bcftools query
- loads into Postgres `variants` table via psql \copy

Designed for offline/local use in docker-compose (`auto_loader` service).
"""

import json
import os
import subprocess
import time
from pathlib import Path

import boto3
from botocore.client import Config


BUCKET = os.getenv("MINIO_DEFAULT_BUCKET", "omics-raw")
PREFIX = os.getenv("MINIO_WATCH_PREFIX", "")
STATE_FILE = Path("/tmp/processed_vcfs.json")
POLL_INTERVAL = int(os.getenv("MINIO_POLL_SECONDS", "60"))

PG_HOST = os.getenv("POSTGRES_HOST", "postgres")
PG_PORT = os.getenv("POSTGRES_PORT", "5432")
PG_USER = os.getenv("POSTGRES_USER", "postgres")
PG_PASSWORD = os.getenv("POSTGRES_PASSWORD", "postgres")
PG_DB = os.getenv("POSTGRES_DB", "omics")
PG_TABLE = os.getenv("POSTGRES_TABLE", "variants")


def get_s3():
    endpoint = os.getenv("MINIO_ENDPOINT", "http://minio:9000")
    access_key = os.getenv("MINIO_ROOT_USER", "minioadmin")
    secret_key = os.getenv("MINIO_ROOT_PASSWORD", "minioadmin123")
    return boto3.client(
        "s3",
        endpoint_url=endpoint,
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        config=Config(signature_version="s3v4"),
        region_name="us-east-1",
    )


def load_state():
    if STATE_FILE.exists():
        return set(json.loads(STATE_FILE.read_text()))
    return set()


def save_state(done):
    STATE_FILE.write_text(json.dumps(sorted(done)))


def list_new_vcfs(s3, done):
    resp = s3.list_objects_v2(Bucket=BUCKET, Prefix=PREFIX)
    contents = resp.get("Contents", [])
    keys = [c["Key"] for c in contents if c["Key"].endswith((".vcf", ".vcf.gz"))]
    return [k for k in keys if k not in done]


def download(s3, key):
    local = Path("/tmp") / key.replace("/", "__")
    local.parent.mkdir(parents=True, exist_ok=True)
    s3.download_file(BUCKET, key, str(local))
    return local


def vcf_to_tsv(vcf_path):
    tsv_path = vcf_path.with_suffix(vcf_path.suffix + ".tsv")
    cmd = [
        "bcftools",
        "query",
        "-H",
        "-f",
        "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\n",
        str(vcf_path),
    ]
    with tsv_path.open("w") as out:
        subprocess.run(cmd, check=True, stdout=out)
    return tsv_path


def load_postgres(tsv_path):
    copy_sql = (
        f"\\copy {PG_TABLE}(chr,pos,ref,alt,qual,filter_status) "
        f"FROM '{tsv_path}' WITH (FORMAT csv, DELIMITER E'\\t', HEADER true)"
    )
    env = os.environ.copy()
    env["PGPASSWORD"] = PG_PASSWORD
    cmd = [
        "psql",
        "--host",
        PG_HOST,
        "--port",
        PG_PORT,
        "--username",
        PG_USER,
        "--dbname",
        PG_DB,
        "-c",
        copy_sql,
    ]
    subprocess.run(cmd, check=True, env=env)


def main():
    s3 = get_s3()
    done = load_state()
    while True:
        try:
            new_keys = list_new_vcfs(s3, done)
            if new_keys:
                print(f"Found {len(new_keys)} new VCF(s): {new_keys}", flush=True)
            for key in new_keys:
                vcf_path = download(s3, key)
                tsv_path = vcf_to_tsv(vcf_path)
                load_postgres(tsv_path)
                done.add(key)
                save_state(done)
                print(f"Ingested {key}", flush=True)
        except Exception as exc:
            print(f"[WARN] Loop error: {exc}", flush=True)
        time.sleep(POLL_INTERVAL)


if __name__ == "__main__":
    main()

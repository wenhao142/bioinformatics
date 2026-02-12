"""
Watch a MinIO bucket for new FASTQ/FASTQ.GZ uploads, run a lightweight
FASTQ -> BAM -> VCF pipeline, upload the VCF back to MinIO, and load
variants into Postgres.

Requirements inside container:
- minimap2, samtools, bcftools, tabix
- Postgres client (psql)
- reference fasta available on local path (REF_FASTA env)
"""

import json
import os
import subprocess
import time
from pathlib import Path
from typing import List, Tuple

import boto3
from botocore.client import Config


BUCKET = os.getenv("MINIO_DEFAULT_BUCKET", "omics-raw")
PREFIX = os.getenv("MINIO_WATCH_PREFIX_FASTQ", "")
STATE_FILE = Path("/tmp/processed_fastq.json")
POLL_INTERVAL = int(os.getenv("MINIO_POLL_SECONDS_FASTQ", "120"))
REF_FASTA = os.getenv("REF_FASTA", "/app/data/GRCh38.fa")

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


def list_fastqs(s3):
    resp = s3.list_objects_v2(Bucket=BUCKET, Prefix=PREFIX)
    contents = resp.get("Contents", [])
    keys = [c["Key"] for c in contents if c["Key"].endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))]
    return keys


def pair_fastqs(keys: List[str]) -> List[Tuple[str, str]]:
    paired = []
    used = set()
    # simple pairing by R1/R2 or _1/_2
    for k in keys:
        if k in used:
            continue
        if "_R1" in k:
            mate = k.replace("_R1", "_R2")
        elif "_1.fastq" in k:
            mate = k.replace("_1.fastq", "_2.fastq")
        elif "_1.fq" in k:
            mate = k.replace("_1.fq", "_2.fq")
        elif "_1.fastq.gz" in k:
            mate = k.replace("_1.fastq.gz", "_2.fastq.gz")
        elif "_1.fq.gz" in k:
            mate = k.replace("_1.fq.gz", "_2.fq.gz")
        else:
            # treat as single-end
            paired.append((k, None))
            used.add(k)
            continue
        if mate in keys:
            paired.append((k, mate))
            used.add(k)
            used.add(mate)
        else:
            # no mate found, treat as single
            paired.append((k, None))
            used.add(k)
    return paired


def download(s3, key):
    local = Path("/tmp") / key.replace("/", "__")
    local.parent.mkdir(parents=True, exist_ok=True)
    s3.download_file(BUCKET, key, str(local))
    return local


def run_pipeline(r1: Path, r2: Path | None, prefix: str):
    bam = Path(f"/tmp/{prefix}.bam")
    vcf = Path(f"/tmp/{prefix}.vcf.gz")

    align_cmd = [
        "minimap2",
        "-t",
        os.getenv("CALLER_THREADS", "4"),
        "-a",
        "-x",
        "sr",
        REF_FASTA,
        str(r1),
    ]
    if r2:
        align_cmd.append(str(r2))

    sort_cmd = ["samtools", "sort", "-o", str(bam)]
    mpileup_cmd = ["bcftools", "mpileup", "-Ou", "-f", REF_FASTA, str(bam)]
    call_cmd = ["bcftools", "call", "-mv", "-Oz", "-o", str(vcf)]
    tabix_cmd = ["tabix", "-f", str(vcf)]

    # Align + sort
    p1 = subprocess.Popen(align_cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p2.communicate()
    subprocess.run(["samtools", "index", str(bam)], check=True)
    # Call variants
    p3 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
    p4 = subprocess.Popen(call_cmd, stdin=p3.stdout)
    p3.stdout.close()
    p4.communicate()
    subprocess.run(tabix_cmd, check=True)
    return vcf


def vcf_to_tsv(vcf_path: Path):
    tsv_path = vcf_path.with_suffix(".tsv")
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


def load_postgres(tsv_path: Path):
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


def upload_vcf(s3, local_vcf: Path, remote_prefix: str):
    key = f"{remote_prefix}{local_vcf.name}"
    s3.upload_file(str(local_vcf), BUCKET, key)
    # upload index
    idx = Path(str(local_vcf) + ".tbi")
    if idx.exists():
        s3.upload_file(str(idx), BUCKET, key + ".tbi")
    return key


def main():
    s3 = get_s3()
    done = load_state()
    while True:
        try:
            keys = list_fastqs(s3)
            pairs = pair_fastqs(keys)
            for r1_key, r2_key in pairs:
                if r1_key in done:
                    continue
                r1_path = download(s3, r1_key)
                r2_path = download(s3, r2_key) if r2_key else None
                prefix = Path(r1_key).stem.replace(".fastq", "").replace(".fq", "")
                vcf_path = run_pipeline(r1_path, r2_path, prefix)
                tsv_path = vcf_to_tsv(vcf_path)
                load_postgres(tsv_path)
                upload_vcf(s3, vcf_path, "results/")
                done.add(r1_key)
                if r2_key:
                    done.add(r2_key)
                save_state(done)
                print(f"Ingested FASTQ -> VCF for {r1_key}", flush=True)
        except Exception as exc:
            print(f"[WARN] Loop error: {exc}", flush=True)
        time.sleep(POLL_INTERVAL)


if __name__ == "__main__":
    main()

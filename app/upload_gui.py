"""
Simple local GUI for uploading genomics files to MinIO (S3-compatible).
- Uses Streamlit; runs inside docker-compose service `uploader`.
- Lets users pick bucket/folder and upload files; lists existing objects for quick verification.
"""

import os
import io
import streamlit as st
import boto3
from botocore.client import Config
from botocore.exceptions import ClientError


def get_s3_client():
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


def ensure_bucket(s3, bucket):
    try:
        s3.head_bucket(Bucket=bucket)
    except ClientError as exc:
        code = int(exc.response.get("ResponseMetadata", {}).get("HTTPStatusCode", 0))
        if code == 404:
            s3.create_bucket(Bucket=bucket)
        else:
            raise


def list_objects(s3, bucket, prefix=""):
    try:
        resp = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    except ClientError:
        return []
    contents = resp.get("Contents", [])
    return [
        {
            "key": obj["Key"],
            "size_MB": round(obj["Size"] / 1024 / 1024, 3),
            "last_modified": obj["LastModified"],
        }
        for obj in contents
    ]


def main():
    st.title("Local Omics Uploader")
    st.caption("Upload your VCF/BCF/FASTQ/CSV to MinIO for downstream Snakemake/analysis.")

    default_bucket = os.getenv("MINIO_DEFAULT_BUCKET", "omics-raw")
    s3 = get_s3_client()

    bucket = st.text_input("Bucket", value=default_bucket)
    prefix = st.text_input("Folder prefix (optional)", value="project1/")

    uploaded = st.file_uploader(
        "Select files", type=None, accept_multiple_files=True, help="Any file type; stored as-is in MinIO."
    )

    if st.button("Upload", type="primary"):
        if not uploaded:
            st.warning("Please choose at least one file.")
        else:
            ensure_bucket(s3, bucket)
            for f in uploaded:
                key = f"{prefix}{f.name}" if prefix else f.name
                data = io.BytesIO(f.getbuffer())
                s3.upload_fileobj(data, bucket, key)
            st.success(f"Uploaded {len(uploaded)} file(s) to s3://{bucket}/{prefix}")

    st.divider()
    st.subheader("Bucket contents")
    objects = []
    if st.button("Refresh list"):
        ensure_bucket(s3, bucket)
        objects = list_objects(s3, bucket, prefix)
    if objects:
        st.dataframe(objects, use_container_width=True)
    else:
        st.write("No objects listed yet. Click Refresh after uploading.")


if __name__ == "__main__":
    main()

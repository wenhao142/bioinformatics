import os
from pathlib import Path
from fastapi.testclient import TestClient

from api.main import app


def get_token(client: TestClient, email: str, password: str) -> str:
    resp = client.post("/auth/login", json={"email": email, "password": password})
    assert resp.status_code == 200
    return resp.json()["access_token"]


def test_upload_uses_local_when_disabled(tmp_path, monkeypatch):
    # force local mode
    monkeypatch.setenv("USE_MINIO", "false")
    monkeypatch.setenv("LOCAL_DATASET_DIR", str(tmp_path))
    client = TestClient(app)

    token = get_token(client, "admin@example.com", "password")
    resp = client.post(
        "/datasets/upload",
        files={"file": ("hello.txt", b"hello-world")},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    rec = resp.json()["dataset"]
    assert rec["filename"] == "hello.txt"
    # verify file written locally
    saved = os.path.join(str(tmp_path), rec["sha256"], "hello.txt")
    assert os.path.exists(saved)


def test_list_datasets_scope(monkeypatch, tmp_path):
    monkeypatch.setenv("USE_MINIO", "false")
    monkeypatch.setenv("LOCAL_DATASET_DIR", str(tmp_path))
    client = TestClient(app)

    admin_token = get_token(client, "admin@example.com", "password")
    viewer_token = get_token(client, "viewer@example.com", "password")

    client.post(
        "/datasets/upload",
        files={"file": ("a.txt", b"aaa")},
        headers={"Authorization": f"Bearer {admin_token}"},
    )
    client.post(
        "/datasets/upload",
        files={"file": ("b.txt", b"bbb")},
        headers={"Authorization": f"Bearer {viewer_token}"},
    )

    # admin sees both
    resp_admin = client.get("/datasets", headers={"Authorization": f"Bearer {admin_token}"})
    assert len(resp_admin.json()["datasets"]) >= 2

    # viewer sees only own
    resp_viewer = client.get("/datasets", headers={"Authorization": f"Bearer {viewer_token}"})
    names = [d["filename"] for d in resp_viewer.json()["datasets"]]
    assert "b.txt" in names
    assert "a.txt" not in names


def test_upload_generated_bio_sample_files(monkeypatch, tmp_path):
    monkeypatch.setenv("USE_MINIO", "false")
    monkeypatch.setenv("LOCAL_DATASET_DIR", str(tmp_path))
    client = TestClient(app)

    token = get_token(client, "admin@example.com", "password")
    headers = {"Authorization": f"Bearer {token}"}

    root = Path(__file__).resolve().parents[2]
    sample_dir = root / "infra" / "sample_data"
    sample_files = [
        "demo_reference.fasta",
        "demo_annotation.gtf",
        "demo_reads_R1.fastq",
        "demo_reads_R2.fastq",
        "demo_samplesheet.tsv",
        "demo_counts.tsv",
        "demo_diff.tsv",
    ]

    uploaded = []
    for name in sample_files:
        path = sample_dir / name
        assert path.exists(), f"Missing sample file: {path}"
        with path.open("rb") as handle:
            resp = client.post(
                "/datasets/upload?project_id=demo-bio",
                files={"file": (name, handle.read())},
                headers=headers,
            )
        assert resp.status_code == 200, resp.text
        body = resp.json()
        assert body["dataset"]["filename"] == name
        uploaded.append(body["dataset"]["filename"])

    assert set(uploaded) == set(sample_files)

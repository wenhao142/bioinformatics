import io

from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient):
    resp = client.post(
        "/auth/login",
        json={"email": "admin@example.com", "password": "password"},
    )
    token = resp.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def ingest_variants(client: TestClient, headers: dict[str, str]):
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n"
        "chr1\t100\t.\tA\tG\t50\tPASS\n"
        "chr1\t130\t.\tA\tT\t40\tPASS\n"
        "chr1\t700\t.\tC\tT\t15\tPASS\n"
    )
    files = {"file": ("run.vcf", io.BytesIO(vcf.encode("utf-8")), "text/vcf")}
    r = client.post("/variants/ingest", headers=headers, files=files)
    assert r.status_code == 200


def upload_omics(client: TestClient, headers: dict[str, str], logfc: float):
    tsv = f"gene\tlogfc\tpval\tadj_pval\nGENE1\t{logfc}\t0.001\t0.005\nGENE2\t0.2\t0.2\t0.4\n"
    files = {"file": ("expr.tsv", io.BytesIO(tsv.encode("utf-8")), "text/tab-separated-values")}
    r = client.post("/omics/expr/upload", headers=headers, files=files)
    assert r.status_code == 200


def test_run_metadata_and_rerun_stability():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)
    upload_omics(client, headers, logfc=1.2)

    first = client.post("/runs/evidence?chr=chr1&start=1&end=1000&top_n=10", headers=headers)
    assert first.status_code == 200
    first_body = first.json()
    first_run = first_body["run"]

    assert first_run["params"] == {"chr": "chr1", "start": 1, "end": 1000, "top_n": 10}
    assert len(first_run["input_hashes"]["variants_region_sha256"]) == 64
    assert len(first_run["input_hashes"]["omics_table_sha256"]) == 64
    assert "python" in first_run["tool_versions"]
    assert first_run["tool_versions"]["scoring"] == "evidence-v1"
    assert first_run["stable_with_previous"] is False
    assert first_run["rerun_of"] is None

    second = client.post("/runs/evidence?chr=chr1&start=1&end=1000&top_n=10", headers=headers)
    assert second.status_code == 200
    second_body = second.json()
    second_run = second_body["run"]

    assert second_run["stable_with_previous"] is True
    assert second_run["rerun_of"] == first_run["run_id"]
    assert second_run["result_hash"] == first_run["result_hash"]

    lookup = client.get(f"/runs/{second_run['run_id']}", headers=headers)
    assert lookup.status_code == 200
    assert lookup.json()["run"]["run_id"] == second_run["run_id"]


def test_rerun_changes_when_inputs_change():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)
    upload_omics(client, headers, logfc=1.2)

    first = client.post("/runs/evidence?chr=chr1&start=1&end=1000", headers=headers)
    assert first.status_code == 200
    first_hash = first.json()["run"]["result_hash"]

    upload_omics(client, headers, logfc=2.4)
    second = client.post("/runs/evidence?chr=chr1&start=1&end=1000", headers=headers)
    assert second.status_code == 200
    second_run = second.json()["run"]

    assert second_run["stable_with_previous"] is False
    assert second_run["rerun_of"] is None
    assert second_run["result_hash"] != first_hash

import io

from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient):
    response = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    token = response.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def ingest_variants(client: TestClient, headers: dict[str, str]):
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n"
        "chr1\t80\t.\tA\tG\t12\tPASS\n"
        "chr1\t110\t.\tC\tT\t20\tPASS\n"
        "chr1\t140\t.\tG\tA\t35\tPASS\n"
        "chr1\t700\t.\tT\tC\t18\tPASS\n"
        "chr1\t760\t.\tA\tG\t28\tPASS\n"
    )
    files = {"file": ("causal.vcf", io.BytesIO(vcf.encode("utf-8")), "text/vcf")}
    response = client.post("/variants/ingest", headers=headers, files=files)
    assert response.status_code == 200


def upload_omics(client: TestClient, headers: dict[str, str]):
    tsv = "gene\tlogfc\tpval\tadj_pval\nGENE1\t1.2\t0.001\t0.005\nGENE2\t0.8\t0.01\t0.03\n"
    files = {"file": ("expr.tsv", io.BytesIO(tsv.encode("utf-8")), "text/tab-separated-values")}
    response = client.post("/omics/expr/upload", headers=headers, files=files)
    assert response.status_code == 200


def test_causal_scoring_run_and_queries():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)
    upload_omics(client, headers)

    run_response = client.post(
        "/causal/score?chr=chr1&start=1&end=1000&top_n=3&ld_window_bp=50000&project_id=demo-project",
        headers=headers,
    )
    assert run_response.status_code == 200
    body = run_response.json()

    run = body["run"]
    result = body["result"]
    assert run["kind"] == "causal_score"
    assert run["project_id"] == "demo-project"
    assert run["counts"]["variant_scores"] == 3
    assert run["counts"]["gene_scores"] >= 1
    assert run["inference_label"] == "inference"
    assert result["meta"]["ld_window_bp"] == 50000
    assert len(result["variant_scores"]) == 3
    assert len(result["gene_scores"]) >= 1
    assert "score" in result["variant_scores"][0]
    assert "score" in result["gene_scores"][0]

    list_response = client.get("/causal/runs", headers=headers)
    assert list_response.status_code == 200
    assert list_response.json()["runs"][0]["run_id"] == run["run_id"]

    run_lookup = client.get(f"/causal/runs/{run['run_id']}", headers=headers)
    assert run_lookup.status_code == 200
    assert run_lookup.json()["run"]["result_hash"] == run["result_hash"]

    result_lookup = client.get(f"/causal/runs/{run['run_id']}/result", headers=headers)
    assert result_lookup.status_code == 200
    assert result_lookup.json()["result"]["meta"]["variants_considered"] >= 5


def test_causal_scoring_accepts_missing_omics():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)

    run_response = client.post("/causal/score?chr=chr1&start=1&end=1000&top_n=2", headers=headers)
    assert run_response.status_code == 200
    body = run_response.json()

    assert body["result"]["meta"]["omics_rows_loaded"] == 0
    assert len(body["result"]["variant_scores"]) == 2

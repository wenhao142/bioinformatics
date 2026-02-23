import io

from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient):
    response = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    token = response.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def seed_inputs(client: TestClient, headers: dict[str, str]):
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n"
        "chr1\t110\t.\tA\tG\t20\tPASS\n"
        "chr1\t140\t.\tC\tT\t30\tPASS\n"
        "chr1\t760\t.\tG\tA\t25\tPASS\n"
    )
    files = {"file": ("report.vcf", io.BytesIO(vcf.encode("utf-8")), "text/vcf")}
    response = client.post("/variants/ingest", headers=headers, files=files)
    assert response.status_code == 200

    tsv = "gene\tlogfc\tpval\tadj_pval\nGENE1\t1.1\t0.001\t0.01\nGENE2\t0.6\t0.03\t0.04\n"
    files = {"file": ("expr.tsv", io.BytesIO(tsv.encode("utf-8")), "text/tab-separated-values")}
    response = client.post("/omics/expr/upload", headers=headers, files=files)
    assert response.status_code == 200


def test_report_markdown_and_html():
    client = TestClient(app)
    headers = auth_header(client)
    seed_inputs(client, headers)

    run_response = client.post("/runs/evidence?chr=chr1&start=1&end=1000&top_n=5&project_id=demo-project", headers=headers)
    assert run_response.status_code == 200
    run_id = run_response.json()["run"]["run_id"]

    markdown = client.get(f"/report/demo-project/{run_id}", headers=headers)
    assert markdown.status_code == 200
    assert markdown.headers["content-type"].startswith("text/markdown")
    assert "## Top Genes" in markdown.text
    assert "## Reproducibility Metadata" in markdown.text

    html = client.get(f"/report/demo-project/{run_id}?format=html", headers=headers)
    assert html.status_code == 200
    assert html.headers["content-type"].startswith("text/html")
    assert "<html>" in html.text
    assert "Analysis Report" in html.text


def test_report_project_mismatch_returns_404():
    client = TestClient(app)
    headers = auth_header(client)
    seed_inputs(client, headers)

    run_response = client.post("/causal/score?chr=chr1&start=1&end=1000&project_id=p1", headers=headers)
    assert run_response.status_code == 200
    run_id = run_response.json()["run"]["run_id"]

    mismatch = client.get(f"/report/p2/{run_id}", headers=headers)
    assert mismatch.status_code == 404

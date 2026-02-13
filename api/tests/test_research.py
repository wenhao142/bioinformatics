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
        "chr1\t700\t.\tC\tT\t15\tPASS\n"
    )
    files = {"file": ("research.vcf", io.BytesIO(vcf.encode("utf-8")), "text/vcf")}
    r = client.post("/variants/ingest", headers=headers, files=files)
    assert r.status_code == 200


def test_offline_template_research_with_input_genes(monkeypatch):
    monkeypatch.setenv("ONLINE_MODE", "false")
    client = TestClient(app)
    headers = auth_header(client)

    payload = {"disease": "Alzheimer disease", "genes": ["APP", "APOE"], "top_n": 2}
    r = client.post("/research/directions", headers=headers, json=payload)
    assert r.status_code == 200
    body = r.json()

    assert body["mode"] == "offline-template"
    assert body["genes_considered"] == ["APP", "APOE"]
    assert len(body["hypotheses"]) == 2
    assert len(body["experiments"]) == 2
    assert body["hypotheses"][0]["inference_label"] == "inference"
    assert body["experiments"][0]["inference_label"] == "inference"


def test_offline_template_research_from_ranked_region():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)

    payload = {"disease": "Alzheimer disease", "chr": "chr1", "start": 1, "end": 1000, "top_n": 2}
    r = client.post("/research/directions", headers=headers, json=payload)
    assert r.status_code == 200
    body = r.json()

    assert body["mode"] == "offline-template"
    assert body["genes_considered"] == ["GENE1", "GENE2"]
    assert len(body["hypotheses"]) == 2
    assert len(body["experiments"]) == 2


def test_research_requires_gene_or_rankable_region():
    client = TestClient(app)
    headers = auth_header(client)

    r = client.post("/research/directions", headers=headers, json={"top_n": 3})
    assert r.status_code == 400
    assert "No genes provided" in r.json()["detail"]

import io

from fastapi.testclient import TestClient

from api import literature
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


def test_research_summary_offline_with_citations(monkeypatch):
    monkeypatch.setenv("ENABLE_LLM_SUMMARY", "false")
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)
    literature._PUBMED_RECORDS.append(
        {
            "gene": "GENE1",
            "pmid": "12345678",
            "title": "Mock AD evidence for GENE1",
            "year": 2024,
            "source": "pubmed",
        }
    )

    payload = {"disease": "Alzheimer disease", "chr": "chr1", "start": 1, "end": 1000, "top_n": 2}
    response = client.post("/research/summary", headers=headers, json=payload)
    assert response.status_code == 200
    body = response.json()

    assert body["mode"] == "offline-template"
    assert body["inference_label"] == "inference"
    assert len(body["citations"]) >= 2
    assert all("id" in row and "source_type" in row for row in body["citations"])
    assert body["llm_mode_requested"] == "auto"
    assert body["llm_model_used"] is None


def test_research_summary_llm_filters_invalid_citations(monkeypatch):
    monkeypatch.setenv("ENABLE_LLM_SUMMARY", "true")
    monkeypatch.setenv("ONLINE_MODE", "true")
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)

    from api import research

    def _fake_llm(_prompt: str, _model: str | None = None):
        return {"summary": "LLM summary", "citations": ["C1", "C999", "C2"]}

    monkeypatch.setattr(research, "_call_cloud_llm", _fake_llm)
    payload = {"disease": "Alzheimer disease", "chr": "chr1", "start": 1, "end": 1000, "top_n": 2}
    response = client.post("/research/summary", headers=headers, json=payload)
    assert response.status_code == 200
    body = response.json()

    assert body["mode"] == "online-llm"
    assert body["summary"] == "LLM summary"
    assert body["citation_ids"] == ["C1", "C2"]
    assert body["llm_mode_requested"] == "auto"


def test_research_summary_forces_offline_even_if_llm_enabled(monkeypatch):
    monkeypatch.setenv("ENABLE_LLM_SUMMARY", "true")
    monkeypatch.setenv("ONLINE_MODE", "true")
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)

    payload = {
        "disease": "Alzheimer disease",
        "chr": "chr1",
        "start": 1,
        "end": 1000,
        "top_n": 2,
        "llm_mode": "offline",
    }
    response = client.post("/research/summary", headers=headers, json=payload)
    assert response.status_code == 200
    body = response.json()

    assert body["mode"] == "offline-template"
    assert body["llm_mode_requested"] == "offline"
    assert body["llm_model_used"] is None


def test_research_summary_uses_model_override(monkeypatch):
    monkeypatch.setenv("ENABLE_LLM_SUMMARY", "true")
    monkeypatch.setenv("ONLINE_MODE", "true")
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)

    from api import research

    calls: list[str | None] = []

    def _fake_llm(_prompt: str, model: str | None = None):
        calls.append(model)
        return {"summary": "Model override summary", "citations": ["C1"]}

    monkeypatch.setattr(research, "_call_cloud_llm", _fake_llm)
    payload = {
        "disease": "Alzheimer disease",
        "chr": "chr1",
        "start": 1,
        "end": 1000,
        "top_n": 2,
        "llm_mode": "auto",
        "llm_model": "gpt-4o-mini",
    }
    response = client.post("/research/summary", headers=headers, json=payload)
    assert response.status_code == 200
    body = response.json()

    assert body["mode"] == "online-llm"
    assert body["llm_model_used"] == "gpt-4o-mini"
    assert calls == ["gpt-4o-mini"]

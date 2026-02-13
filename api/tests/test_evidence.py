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
        "chr1\t760\t.\tG\tA\t10\tPASS\n"
    )
    files = {"file": ("rank.vcf", io.BytesIO(vcf.encode("utf-8")), "text/vcf")}
    r = client.post("/variants/ingest", headers=headers, files=files)
    assert r.status_code == 200
    assert r.json()["ingested"] == 4


def test_rank_evidence_with_feature_breakdown():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)

    tsv = (
        "gene\tlogfc\tpval\tadj_pval\n"
        "GENE1\t1.2\t0.001\t0.005\n"
        "GENE2\t0.2\t0.2\t0.4\n"
    )
    files = {"file": ("expr.tsv", io.BytesIO(tsv.encode("utf-8")), "text/tab-separated-values")}
    upload = client.post("/omics/expr/upload", headers=headers, files=files)
    assert upload.status_code == 200

    r = client.get("/evidence/rank?chr=chr1&start=1&end=1000&top_n=10", headers=headers)
    assert r.status_code == 200
    body = r.json()

    assert body["meta"]["variants_considered"] == 4
    assert len(body["ranked_genes"]) >= 2
    assert body["ranked_genes"][0]["gene"] == "GENE1"
    assert body["ranked_genes"][0]["feature_breakdown"]["variant_count"] == 2
    assert "components" in body["ranked_genes"][0]["feature_breakdown"]

    assert len(body["ranked_loci"]) >= 2
    assert body["ranked_loci"][0]["gene"] == "GENE1"
    assert body["ranked_loci"][0]["feature_breakdown"]["omics"]["adj_pval"] == 0.005


def test_rank_evidence_without_omics_rows():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_variants(client, headers)

    r = client.get("/evidence/rank?chr=chr1&start=1&end=1000", headers=headers)
    assert r.status_code == 200
    body = r.json()

    assert body["meta"]["omics_rows_loaded"] == 0
    assert len(body["ranked_genes"]) >= 2
    top = body["ranked_genes"][0]["feature_breakdown"]["omics"]
    assert top["logfc"] is None
    assert top["effective_pval"] is None

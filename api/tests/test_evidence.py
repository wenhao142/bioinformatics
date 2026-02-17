import io

from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient):
    resp = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    token = resp.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def test_aggregate_with_variants_and_expr():
    client = TestClient(app)
    headers = auth_header(client)

    # upload expression
    tsv = "gene\tlogfc\tpval\nGENE1\t2.0\t0.001\nGENE2\t-1.0\t0.05\n"
    files = {"file": ("expr.tsv", io.BytesIO(tsv.encode()), "text/tab-separated-values")}
    r = client.post("/omics/expr/upload", headers=headers, files=files)
    assert r.status_code == 200

    # ingest variants that map to GENE1/GENE2 intervals
    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\nchr1\t50\t.\tA\tG\t10\tPASS\nchr1\t700\t.\tC\tT\t20\tPASS\n"
    files = {"file": ("simple.vcf", io.BytesIO(vcf.encode()), "text/vcf")}
    r = client.post("/variants/ingest", headers=headers, files=files)
    assert r.status_code == 200

    r = client.get("/evidence/aggregate", headers=headers)
    assert r.status_code == 200
    genes = r.json()["genes"]
    assert genes[0]["gene"] in {"GENE1", "GENE2"}
    # should contain scores and counts
    assert all("score" in g and "variant_count" in g for g in genes)

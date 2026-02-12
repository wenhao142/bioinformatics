from fastapi.testclient import TestClient

from api.main import app


def get_token(client: TestClient, email: str, password: str) -> str:
    resp = client.post("/auth/login", json={"email": email, "password": password})
    assert resp.status_code == 200
    return resp.json()["access_token"]


VCF_CONTENT = b"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t50\tPASS\t.
chr1\t150\t.\tT\tC,A\t20\tq10\t.
"""


def test_ingest_and_query(monkeypatch):
    # use memory store
    monkeypatch.setenv("VARIANTS_USE_DB", "false")
    client = TestClient(app)
    token = get_token(client, "admin@example.com", "password")

    resp = client.post(
        "/variants/ingest",
        files={"file": ("small.vcf", VCF_CONTENT, "text/plain")},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    assert resp.json()["ingested"] == 3  # two lines, three alts

    query = client.get(
        "/variants",
        params={"chr": "chr1", "start": 90, "end": 200},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert query.status_code == 200
    assert len(query.json()["variants"]) == 3

    stats = client.get("/variants/stats", headers={"Authorization": f"Bearer {token}"})
    assert stats.status_code == 200
    assert stats.json()["total_variants"] == 3

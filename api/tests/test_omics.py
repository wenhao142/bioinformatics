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


def test_upload_and_stats():
    client = TestClient(app)
    hdrs = auth_header(client)

    tsv = "gene\tlogfc\tpval\tadj_pval\nA\t1.0\t0.01\t0.02\nB\t-0.5\t0.2\t0.2\n"
    files = {"file": ("expr.tsv", io.BytesIO(tsv.encode("utf-8")), "text/tab-separated-values")}
    r = client.post("/omics/expr/upload", headers=hdrs, files=files)
    assert r.status_code == 200
    assert r.json()["ingested"] == 2

    r = client.get("/omics/expr", headers=hdrs)
    assert r.status_code == 200
    assert len(r.json()["rows"]) == 2

    r = client.get("/omics/expr/stats", headers=hdrs)
    assert r.status_code == 200
    body = r.json()
    assert body["count"] == 2
    assert body["sig_genes"] == 1

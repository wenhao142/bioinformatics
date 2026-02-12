import io
import os

# Force memory store for tests to keep isolation
os.environ["VARIANTS_USE_DB"] = "false"

from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient):
    resp = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    token = resp.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def ingest_simple(client, headers):
    # reset in-memory store if present to keep tests isolated
    try:
        import api.api.variants as vmod

        if hasattr(vmod.STORE, "rows"):
            vmod.STORE.rows.clear()
        elif hasattr(vmod.STORE, "conn"):
            with vmod.STORE.conn.cursor() as cur:
                cur.execute("DELETE FROM variants")
    except Exception:
        pass
    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\nchr1\t50\t.\tA\tG\t10\tPASS\nchr1\t800\t.\tC\tT\t20\tPASS\n"
    files = {"file": ("simple.vcf", io.BytesIO(vcf.encode()), "text/vcf")}
    r = client.post("/variants/ingest", headers=headers, files=files)
    assert r.status_code == 200


def test_nearest_gene_endpoint():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_simple(client, headers)

    r = client.get("/variants/nearest?chr=chr1&pos=60", headers=headers)
    assert r.status_code == 200
    assert r.json()["gene"] == "GENE1"

    r = client.get("/variants/nearest?chr=chr1&pos=850", headers=headers)
    assert r.status_code == 200
    assert r.json()["gene"] == "GENE2"


def test_list_with_annotation():
    client = TestClient(app)
    headers = auth_header(client)
    ingest_simple(client, headers)

    r = client.get("/variants?chr=chr1&start=1&end=1000&annotate=true", headers=headers)
    body = r.json()
    assert len(body["variants"]) >= 2
    genes = {v["gene"] for v in body["variants"]}
    assert {"GENE1", "GENE2"}.issubset(genes)

import io

from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient):
    resp = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    token = resp.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def test_bed_export_memory_store():
    client = TestClient(app)
    headers = auth_header(client)

    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\nchr1\t50\t.\tA\tG\t10\tPASS\nchr1\t800\t.\tC\tT\t20\tPASS\n"
    files = {"file": ("simple.vcf", io.BytesIO(vcf.encode()), "text/vcf")}
    r = client.post("/variants/ingest", headers=headers, files=files)
    assert r.status_code == 200

    r = client.get("/variants/bed?chr=chr1&start=1&end=1000", headers=headers)
    assert r.status_code == 200
    lines = [ln for ln in r.text.splitlines() if ln.strip()]
    assert len(lines) == 2
    assert lines[0].startswith("chr1\t49\t50\tA>G")


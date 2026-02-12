from fastapi.testclient import TestClient

from api.main import app


def get_token(client: TestClient, email: str, password: str) -> str:
    resp = client.post("/auth/login", json={"email": email, "password": password})
    assert resp.status_code == 200
    return resp.json()["access_token"]


def test_audit_records_upload_and_export():
    client = TestClient(app)
    admin_token = get_token(client, "admin@example.com", "password")

    # create project
    client.post(
        "/projects",
        json={"project_id": "p-audit"},
        headers={"Authorization": f"Bearer {admin_token}"},
    )
    # run upload + export
    client.post(
        "/projects/p-audit/upload",
        headers={"Authorization": f"Bearer {admin_token}"},
    )
    client.post(
        "/projects/p-audit/export",
        headers={"Authorization": f"Bearer {admin_token}"},
    )

    resp = client.get("/audit", headers={"Authorization": f"Bearer {admin_token}"})
    assert resp.status_code == 200
    actions = [e["action"] for e in resp.json()["events"] if e["project_id"] == "p-audit"]
    assert "upload" in actions
    assert "export" in actions

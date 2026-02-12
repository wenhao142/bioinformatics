from fastapi.testclient import TestClient

from api.main import app


def get_token(client: TestClient, email: str, password: str) -> str:
    resp = client.post("/auth/login", json={"email": email, "password": password})
    assert resp.status_code == 200, resp.text
    return resp.json()["access_token"]


def test_viewer_cannot_upload(monkeypatch):
    client = TestClient(app)
    # Seed project with admin and viewer
    token_admin = get_token(client, "admin@example.com", "password")

    # Create project
    client.post(
        "/projects",
        json={"project_id": "p1"},
        headers={"Authorization": f"Bearer {token_admin}"},
    )
    # Add viewer
    client.post(
        "/projects/p1/members",
        json={"email": "viewer@example.com", "role": "viewer"},
        headers={"Authorization": f"Bearer {token_admin}"},
    )
    viewer_token = get_token(client, "viewer@example.com", "password")

    # Viewer can view
    resp_view = client.get(
        "/projects/p1/view",
        headers={"Authorization": f"Bearer {viewer_token}"},
    )
    assert resp_view.status_code == 200

    # Viewer cannot upload
    resp_upload = client.post(
        "/projects/p1/upload",
        headers={"Authorization": f"Bearer {viewer_token}"},
    )
    assert resp_upload.status_code == 403


def test_analyst_can_upload():
    client = TestClient(app)
    token_admin = get_token(client, "admin@example.com", "password")
    client.post(
        "/projects",
        json={"project_id": "p2"},
        headers={"Authorization": f"Bearer {token_admin}"},
    )
    client.post(
        "/projects/p2/members",
        json={"email": "analyst@example.com", "role": "analyst"},
        headers={"Authorization": f"Bearer {token_admin}"},
    )
    analyst_token = get_token(client, "analyst@example.com", "password")
    resp = client.post(
        "/projects/p2/upload",
        headers={"Authorization": f"Bearer {analyst_token}"},
    )
    assert resp.status_code == 200

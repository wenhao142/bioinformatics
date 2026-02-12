from fastapi.testclient import TestClient

from api.main import app


def test_login_success_and_me():
    client = TestClient(app)
    resp = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    assert resp.status_code == 200
    token = resp.json()["access_token"]

    me = client.get("/auth/me", headers={"Authorization": f"Bearer {token}"})
    assert me.status_code == 200
    assert me.json()["email"] == "admin@example.com"


def test_login_failure():
    client = TestClient(app)
    resp = client.post("/auth/login", json={"email": "admin@example.com", "password": "wrong"})
    assert resp.status_code == 401

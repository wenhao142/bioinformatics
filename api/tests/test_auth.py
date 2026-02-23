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


def test_admin_can_register_user_and_new_user_can_login():
    client = TestClient(app)
    admin_login = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    assert admin_login.status_code == 200
    admin_token = admin_login.json()["access_token"]

    register = client.post(
        "/auth/register",
        headers={"Authorization": f"Bearer {admin_token}"},
        json={"email": "new.user@example.com", "password": "newpassword1"},
    )
    assert register.status_code == 201
    assert register.json()["role"] == "admin"

    login = client.post("/auth/login", json={"email": "new.user@example.com", "password": "newpassword1"})
    assert login.status_code == 200


def test_non_admin_cannot_register_user():
    client = TestClient(app)
    analyst_login = client.post("/auth/login", json={"email": "analyst@example.com", "password": "password"})
    assert analyst_login.status_code == 200
    analyst_token = analyst_login.json()["access_token"]

    register = client.post(
        "/auth/register",
        headers={"Authorization": f"Bearer {analyst_token}"},
        json={"email": "blocked@example.com", "password": "newpassword1", "role": "viewer"},
    )
    assert register.status_code == 403

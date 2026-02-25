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
        json={"email": "new.user@example.com", "password": "newpassword1", "role": "analyst"},
    )
    assert register.status_code == 201
    assert register.json()["role"] == "analyst"

    login = client.post("/auth/login", json={"email": "new.user@example.com", "password": "newpassword1"})
    assert login.status_code == 200
    me = client.get("/auth/me", headers={"Authorization": f"Bearer {login.json()['access_token']}"})
    assert me.status_code == 200
    assert me.json()["role"] == "analyst"


def test_register_defaults_to_viewer_role():
    client = TestClient(app)
    admin_login = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    assert admin_login.status_code == 200
    admin_token = admin_login.json()["access_token"]

    register = client.post(
        "/auth/register",
        headers={"Authorization": f"Bearer {admin_token}"},
        json={"email": "viewer.default@example.com", "password": "newpassword1"},
    )
    assert register.status_code == 201
    assert register.json()["role"] == "viewer"


def test_register_rejects_invalid_role():
    client = TestClient(app)
    admin_login = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    assert admin_login.status_code == 200
    admin_token = admin_login.json()["access_token"]

    register = client.post(
        "/auth/register",
        headers={"Authorization": f"Bearer {admin_token}"},
        json={"email": "bad.role@example.com", "password": "newpassword1", "role": "owner"},
    )
    assert register.status_code == 400


def test_admin_user_list_does_not_expose_password_hash():
    client = TestClient(app)
    admin_login = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    assert admin_login.status_code == 200
    admin_token = admin_login.json()["access_token"]

    resp = client.get("/auth/users", headers={"Authorization": f"Bearer {admin_token}"})
    assert resp.status_code == 200
    users = resp.json()["users"]
    assert len(users) >= 1
    for user in users:
        assert set(user.keys()) == {"email", "role"}
        assert "password_hash" not in user


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

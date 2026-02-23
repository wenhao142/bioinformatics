from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient, email: str, password: str = "password"):
    response = client.post("/auth/login", json={"email": email, "password": password})
    assert response.status_code == 200
    token = response.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def valid_manifest(plugin_id: str = "baseline-causal"):
    return {
        "plugin_id": plugin_id,
        "name": "Baseline Causal Plugin",
        "version": "1.0.0",
        "image": "ghcr.io/example/causal-plugin:1.0.0",
        "description": "Scores variants and genes from region-level evidence.",
        "input_schema": {
            "type": "object",
            "properties": {
                "project_id": {"type": "string"},
                "chr": {"type": "string"},
                "start": {"type": "integer"},
                "end": {"type": "integer"},
            },
            "required": ["project_id", "chr", "start", "end"],
            "additionalProperties": False,
        },
        "output_schema": {
            "type": "object",
            "properties": {
                "run_id": {"type": "string"},
                "gene_scores": {"type": "array"},
                "variant_scores": {"type": "array"},
            },
            "required": ["run_id", "gene_scores", "variant_scores"],
        },
        "resources": {"cpu_millicores": 1000, "memory_mb": 2048, "gpu": False},
        "enabled": True,
        "tags": ["baseline", "offline"],
    }


def test_register_plugin_manifest_success():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    response = client.post("/plugins/register", headers=headers, json=valid_manifest())
    assert response.status_code == 201
    body = response.json()
    assert body["plugin"]["plugin_id"] == "baseline-causal"
    assert body["plugin"]["registered_by"] == "admin@example.com"

    listed = client.get("/plugins", headers=headers)
    assert listed.status_code == 200
    assert len(listed.json()["plugins"]) == 1


def test_register_plugin_schema_validation_error():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = valid_manifest("bad-schema")
    manifest["input_schema"]["required"] = ["missing_key"]
    response = client.post("/plugins/register", headers=headers, json=manifest)
    assert response.status_code == 400
    assert "required contains unknown fields" in response.json()["detail"]


def test_register_plugin_admin_only_and_duplicate():
    client = TestClient(app)
    admin_headers = auth_header(client, "admin@example.com")
    viewer_headers = auth_header(client, "viewer@example.com")

    forbidden = client.post("/plugins/register", headers=viewer_headers, json=valid_manifest("v1"))
    assert forbidden.status_code == 403

    first = client.post("/plugins/register", headers=admin_headers, json=valid_manifest("dup-plugin"))
    assert first.status_code == 201
    second = client.post("/plugins/register", headers=admin_headers, json=valid_manifest("dup-plugin"))
    assert second.status_code == 409

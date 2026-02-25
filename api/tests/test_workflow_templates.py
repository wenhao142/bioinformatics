import json
from pathlib import Path

from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient):
    response = client.post("/auth/login", json={"email": "admin@example.com", "password": "password"})
    assert response.status_code == 200
    token = response.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def _template_paths() -> list[Path]:
    root = Path(__file__).resolve().parents[2]
    templates_dir = root / "workflow" / "templates"
    return sorted(templates_dir.glob("*_template.json"))


def test_all_workflow_templates_can_be_imported_and_validated():
    client = TestClient(app)
    headers = auth_header(client)

    paths = _template_paths()
    assert len(paths) >= 5
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        response = client.post("/workflows/import", headers=headers, json=payload)
        assert response.status_code == 201, f"{path.name}: {response.text}"
        assert response.json()["report"]["acyclic"] is True

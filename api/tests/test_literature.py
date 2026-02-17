from fastapi.testclient import TestClient

from api import literature
from api.main import app


def auth_header(client: TestClient):
    resp = client.post(
        "/auth/login",
        json={"email": "admin@example.com", "password": "password"},
    )
    token = resp.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def test_pubmed_disabled_in_offline_mode(monkeypatch):
    monkeypatch.delenv("ONLINE_MODE", raising=False)
    monkeypatch.delenv("ENABLE_PUBMED", raising=False)

    client = TestClient(app)
    headers = auth_header(client)
    r = client.get("/literature/pubmed/fetch?genes=APP", headers=headers)
    assert r.status_code == 503
    assert "offline mode" in r.json()["detail"]


def test_pubmed_fetch_and_store(monkeypatch):
    monkeypatch.setenv("ONLINE_MODE", "true")

    def fake_http_get_json(url: str):
        if "esearch.fcgi" in url:
            return {"esearchresult": {"idlist": ["12345", "67890"]}}
        if "esummary.fcgi" in url:
            return {
                "result": {
                    "uids": ["12345", "67890"],
                    "12345": {"uid": "12345", "title": "APP in Alzheimer disease", "pubdate": "2024 Jan"},
                    "67890": {"uid": "67890", "title": "APOE and dementia", "pubdate": "2018"},
                }
            }
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(literature, "_http_get_json", fake_http_get_json)

    client = TestClient(app)
    headers = auth_header(client)

    fetch = client.get(
        "/literature/pubmed/fetch?genes=APP&genes=APOE&disease=Alzheimer%20disease&max_per_gene=2",
        headers=headers,
    )
    assert fetch.status_code == 200
    body = fetch.json()
    assert body["mode"] == "online"
    assert len(body["records"]) == 4
    assert body["errors"] == []

    first = body["records"][0]
    assert first["pmid"] == "12345"
    assert first["title"] == "APP in Alzheimer disease"
    assert first["year"] == 2024

    listed = client.get("/literature/pubmed?gene=APOE", headers=headers)
    assert listed.status_code == 200
    listed_body = listed.json()
    assert listed_body["total"] == 2
    assert listed_body["records"][0]["gene"] == "APOE"

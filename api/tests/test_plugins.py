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


def runner_manifest(plugin_id: str = "baseline-runner"):
    manifest = valid_manifest(plugin_id)
    manifest["image"] = "builtin/baseline-evidence:1.0.0"
    manifest["output_schema"] = {
        "type": "object",
        "properties": {
            "region": {"type": "object"},
            "ranked_genes": {"type": "array"},
            "ranked_loci": {"type": "array"},
            "meta": {"type": "object"},
        },
        "required": ["ranked_genes", "ranked_loci"],
    }
    return manifest


def ingest_variants(client: TestClient, headers: dict[str, str]):
    content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n"
        "chr1\t80\t.\tA\tG\t12\tPASS\n"
        "chr1\t110\t.\tC\tT\t20\tPASS\n"
        "chr1\t140\t.\tG\tA\t35\tPASS\n"
        "chr1\t700\t.\tT\tC\t18\tPASS\n"
    )
    files = {"file": ("plugin.vcf", content, "text/vcf")}
    response = client.post("/variants/ingest", headers=headers, files=files)
    assert response.status_code == 200


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


def test_baseline_plugin_run_persists_ranked_outputs():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    ingest_variants(client, headers)

    register = client.post("/plugins/register", headers=headers, json=runner_manifest("baseline-runner"))
    assert register.status_code == 201

    run = client.post(
        "/plugins/baseline-runner/run",
        headers=headers,
        json={"project_id": "demo-project", "chr": "chr1", "start": 1, "end": 1000, "top_n": 3},
    )
    assert run.status_code == 200
    body = run.json()
    assert body["run"]["engine"] == "builtin"
    assert body["run"]["counts"]["ranked_genes"] >= 1
    assert len(body["result"]["ranked_genes"]) >= 1
    run_id = body["run"]["run_id"]

    listed = client.get("/plugins/baseline-runner/runs", headers=headers)
    assert listed.status_code == 200
    assert listed.json()["runs"][0]["run_id"] == run_id

    result = client.get(f"/plugins/baseline-runner/runs/{run_id}/result", headers=headers)
    assert result.status_code == 200
    assert "ranked_loci" in result.json()["result"]

    disabled = client.patch("/plugins/baseline-runner/enabled", headers=headers, json={"enabled": False})
    assert disabled.status_code == 200
    assert disabled.json()["enabled"] is False

    blocked = client.post(
        "/plugins/baseline-runner/run",
        headers=headers,
        json={"project_id": "demo-project", "chr": "chr1", "start": 1, "end": 1000, "top_n": 2},
    )
    assert blocked.status_code == 400

    enabled = client.patch("/plugins/baseline-runner/enabled", headers=headers, json={"enabled": True})
    assert enabled.status_code == 200
    assert enabled.json()["enabled"] is True


def test_docker_plugin_runner_path(monkeypatch):
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = valid_manifest("docker-plugin")
    manifest["tags"] = ["docker"]
    response = client.post("/plugins/register", headers=headers, json=manifest)
    assert response.status_code == 201

    monkeypatch.setenv("PLUGIN_RUNNER_ALLOW_DOCKER", "true")
    from api import plugins

    class _Result:
        returncode = 0
        stdout = '{"run_id":"x","gene_scores":[],"variant_scores":[]}'
        stderr = ""

    def _fake_run(*_args, **_kwargs):
        return _Result()

    monkeypatch.setattr(plugins.subprocess, "run", _fake_run)
    run = client.post(
        "/plugins/docker-plugin/run",
        headers=headers,
        json={"project_id": "demo-project", "chr": "chr1", "start": 1, "end": 1000, "top_n": 2},
    )
    assert run.status_code == 200
    assert run.json()["run"]["engine"] == "docker"


def test_compare_plugin_runs():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    ingest_variants(client, headers)
    register = client.post("/plugins/register", headers=headers, json=runner_manifest("compare-runner"))
    assert register.status_code == 201

    first = client.post(
        "/plugins/compare-runner/run",
        headers=headers,
        json={"project_id": "demo-project", "chr": "chr1", "start": 1, "end": 1000, "top_n": 1},
    )
    assert first.status_code == 200
    first_run_id = first.json()["run"]["run_id"]

    second = client.post(
        "/plugins/compare-runner/run",
        headers=headers,
        json={"project_id": "demo-project", "chr": "chr1", "start": 1, "end": 1000, "top_n": 3},
    )
    assert second.status_code == 200
    second_run_id = second.json()["run"]["run_id"]

    compare = client.post(
        "/plugins/compare-runner/compare",
        headers=headers,
        json={"left_run_id": first_run_id, "right_run_id": second_run_id},
    )
    assert compare.status_code == 200
    body = compare.json()
    assert body["left_run_id"] == first_run_id
    assert body["right_run_id"] == second_run_id
    assert "overlap_genes" in body
    assert "score_deltas" in body

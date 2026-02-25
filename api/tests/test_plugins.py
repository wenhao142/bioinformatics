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


def adapter_manifest(plugin_id: str = "vcf-adapter", version: str = "1.0.0"):
    manifest = valid_manifest(plugin_id)
    manifest["version"] = version
    manifest["tags"] = ["adapter"]
    manifest["adapter"] = {
        "from_type": "variants.vcf.gz",
        "to_type": "variants.vcf.gz",
        "strategy": "normalize",
    }
    return manifest


def wrapper_manifest(plugin_id: str, wrapper_tag: str):
    manifest = valid_manifest(plugin_id)
    manifest["image"] = f"builtin/{wrapper_tag}:1.0.0"
    manifest["tags"] = [wrapper_tag]
    manifest["output_schema"] = {
        "type": "object",
        "properties": {
            "status": {"type": "string"},
            "workflow_engine": {"type": "string"},
            "command_preview": {"type": "object"},
            "ranked_genes": {"type": "array"},
            "ranked_loci": {"type": "array"},
        },
        "required": ["status", "workflow_engine", "command_preview", "ranked_genes", "ranked_loci"],
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


def test_list_canonical_types():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    response = client.get("/plugins/canonical-types", headers=headers)
    assert response.status_code == 200
    body = response.json()
    assert body["registry_version"] == "1.0.0"
    type_ids = {row["id"] for row in body["types"]}
    assert "reads.fastq.gz" in type_ids
    assert "variants.vcf.gz" in type_ids


def test_method_registry_supports_versions_and_latest_resolution():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")

    first = valid_manifest("versioned-plugin")
    first["version"] = "1.0.0"
    second = valid_manifest("versioned-plugin")
    second["version"] = "1.1.0"

    r1 = client.post("/plugins/register", headers=headers, json=first)
    assert r1.status_code == 201
    assert r1.json()["plugin"]["version"] == "1.0.0"
    assert r1.json()["plugin"]["is_latest"] is True

    r2 = client.post("/plugins/register", headers=headers, json=second)
    assert r2.status_code == 201
    assert r2.json()["plugin"]["version"] == "1.1.0"
    assert r2.json()["plugin"]["is_latest"] is True

    latest = client.get("/plugins/versioned-plugin", headers=headers)
    assert latest.status_code == 200
    assert latest.json()["plugin"]["version"] == "1.1.0"
    assert latest.json()["plugin"]["is_latest"] is True

    old_version = client.get("/plugins/versioned-plugin?version=1.0.0", headers=headers)
    assert old_version.status_code == 200
    assert old_version.json()["plugin"]["version"] == "1.0.0"
    assert old_version.json()["plugin"]["is_latest"] is False

    latest_only = client.get("/plugins", headers=headers)
    assert latest_only.status_code == 200
    assert len(latest_only.json()["plugins"]) == 1
    assert latest_only.json()["plugins"][0]["version"] == "1.1.0"

    all_versions = client.get("/plugins?all_versions=true", headers=headers)
    assert all_versions.status_code == 200
    versions = [row["version"] for row in all_versions.json()["plugins"]]
    assert versions == ["1.1.0", "1.0.0"]


def test_register_plugin_schema_validation_error():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = valid_manifest("bad-schema")
    manifest["input_schema"]["required"] = ["missing_key"]
    response = client.post("/plugins/register", headers=headers, json=manifest)
    assert response.status_code == 400
    assert "required contains unknown fields" in response.json()["detail"]


def test_register_plugin_manifest_rejects_unknown_top_level_field():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = valid_manifest("bad-top-field")
    manifest["unknown_field"] = {"x": 1}
    response = client.post("/plugins/register", headers=headers, json=manifest)
    assert response.status_code == 400
    assert "Additional properties are not allowed" in response.json()["detail"]


def test_register_plugin_manifest_rejects_invalid_resources_shape():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = valid_manifest("bad-resources")
    manifest["resources"] = {"cpu_millicores": 50, "memory_mb": 2048, "gpu": False}
    response = client.post("/plugins/register", headers=headers, json=manifest)
    assert response.status_code == 400
    assert "resources.cpu_millicores" in response.json()["detail"]


def test_register_plugin_manifest_rejects_unknown_canonical_type():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = valid_manifest("bad-canonical")
    manifest["input_schema"]["properties"]["chr"]["canonical_type"] = "reads.fastq"
    response = client.post("/plugins/register", headers=headers, json=manifest)
    assert response.status_code == 400
    assert "canonical_type is unknown" in response.json()["detail"]


def test_register_plugin_manifest_accepts_known_canonical_type():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = valid_manifest("good-canonical")
    manifest["input_schema"]["properties"]["chr"]["canonical_type"] = "variants.vcf.gz"
    response = client.post("/plugins/register", headers=headers, json=manifest)
    assert response.status_code == 201


def test_register_plugin_rejects_latest_image_tag():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = valid_manifest("bad-image-tag")
    manifest["image"] = "ghcr.io/example/causal-plugin:latest"
    response = client.post("/plugins/register", headers=headers, json=manifest)
    assert response.status_code == 400
    assert "latest" in response.json()["detail"]


def test_adapter_plugin_requires_spec_and_is_listed():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")

    missing = valid_manifest("missing-adapter-spec")
    missing["tags"] = ["adapter"]
    missing_resp = client.post("/plugins/register", headers=headers, json=missing)
    assert missing_resp.status_code == 400
    assert "adapter plugin requires adapter specification" in missing_resp.json()["detail"]

    good_resp = client.post("/plugins/register", headers=headers, json=adapter_manifest("vcf-adapter"))
    assert good_resp.status_code == 201
    adapters = client.get("/plugins/adapters", headers=headers)
    assert adapters.status_code == 200
    assert any(row["plugin_id"] == "vcf-adapter" for row in adapters.json()["adapters"])


def test_plugin_search_api_filters():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    first = client.post("/plugins/register", headers=headers, json=adapter_manifest("search-adapter"))
    assert first.status_code == 201
    second = client.post("/plugins/register", headers=headers, json=valid_manifest("search-baseline"))
    assert second.status_code == 201

    by_tag = client.get("/plugins/search?tag=adapter", headers=headers)
    assert by_tag.status_code == 200
    ids = [row["plugin_id"] for row in by_tag.json()["plugins"]]
    assert "search-adapter" in ids
    assert "search-baseline" not in ids

    by_query = client.get("/plugins/search?q=baseline", headers=headers)
    assert by_query.status_code == 200
    ids_q = [row["plugin_id"] for row in by_query.json()["plugins"]]
    assert "search-baseline" in ids_q


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
    assert body["run"]["plugin_version"] == "1.0.0"
    assert body["run"]["params_hash"]
    assert body["run"]["tool_versions"]["plugin_runner"]
    assert body["run"]["execution"]["duration_ms"] >= 0
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
        def __init__(self, returncode: int = 0, stdout: str = "", stderr: str = ""):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr

    def _fake_run(args, *_other, **_kwargs):
        cmd = " ".join(args)
        if "image inspect" in cmd:
            return _Result(returncode=0, stdout="ghcr.io/example/causal-plugin@sha256:abc")
        return _Result(returncode=0, stdout='{"run_id":"x","gene_scores":[],"variant_scores":[]}', stderr="")

    monkeypatch.setattr(plugins.subprocess, "run", _fake_run)
    run = client.post(
        "/plugins/docker-plugin/run",
        headers=headers,
        json={"project_id": "demo-project", "chr": "chr1", "start": 1, "end": 1000, "top_n": 2},
    )
    assert run.status_code == 200
    assert run.json()["run"]["engine"] == "docker"
    assert run.json()["run"]["container_image_digest"] == "ghcr.io/example/causal-plugin@sha256:abc"


def test_deprecated_plugin_cannot_run():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = runner_manifest("deprecated-runner")
    manifest["deprecated"] = True
    manifest["deprecation_message"] = "Use baseline-runner-v2"
    register = client.post("/plugins/register", headers=headers, json=manifest)
    assert register.status_code == 201
    run = client.post(
        "/plugins/deprecated-runner/run",
        headers=headers,
        json={"project_id": "demo-project", "chr": "chr1", "start": 1, "end": 1000, "top_n": 2},
    )
    assert run.status_code == 400
    assert "deprecated" in run.json()["detail"].lower()


def test_snakemake_wrapper_plugin_run():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = wrapper_manifest("snakemake-wrapper-plugin", "snakemake-wrapper")
    register = client.post("/plugins/register", headers=headers, json=manifest)
    assert register.status_code == 201
    run = client.post(
        "/plugins/snakemake-wrapper-plugin/run",
        headers=headers,
        json={
            "project_id": "demo-project",
            "chr": "chr1",
            "start": 1,
            "end": 1000,
            "top_n": 2,
            "parameters": {"workflow_file": "workflow/Snakefile", "dry_run": True},
        },
    )
    assert run.status_code == 200
    assert run.json()["run"]["engine"] == "snakemake-wrapper"
    assert run.json()["result"]["workflow_engine"] == "snakemake"


def test_nextflow_wrapper_plugin_run():
    client = TestClient(app)
    headers = auth_header(client, "admin@example.com")
    manifest = wrapper_manifest("nextflow-wrapper-plugin", "nextflow-wrapper")
    register = client.post("/plugins/register", headers=headers, json=manifest)
    assert register.status_code == 201
    run = client.post(
        "/plugins/nextflow-wrapper-plugin/run",
        headers=headers,
        json={
            "project_id": "demo-project",
            "chr": "chr1",
            "start": 1,
            "end": 1000,
            "top_n": 2,
            "parameters": {"workflow_file": "main.nf", "dry_run": True},
        },
    )
    assert run.status_code == 200
    assert run.json()["run"]["engine"] == "nextflow-wrapper"
    assert run.json()["result"]["workflow_engine"] == "nextflow"


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

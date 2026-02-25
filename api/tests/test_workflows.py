from fastapi.testclient import TestClient

from api.main import app


def auth_header(client: TestClient, email: str = "admin@example.com", password: str = "password"):
    response = client.post("/auth/login", json={"email": email, "password": password})
    assert response.status_code == 200
    token = response.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def valid_workflow(workflow_id: str = "wf-branch"):
    return {
        "workflow_id": workflow_id,
        "nodes": [
            {
                "node_id": "n1",
                "plugin_id": "source",
                "output_types": {"vcf": "variants.vcf.gz"},
                "parameters": {"mode": "source"},
            },
            {
                "node_id": "n2",
                "plugin_id": "consumer-a",
                "input_types": {"input_vcf": "variants.vcf.gz"},
                "output_types": {"rep": "report.html"},
            },
            {
                "node_id": "n3",
                "plugin_id": "consumer-b",
                "input_types": {"input_vcf": "variants.vcf.gz"},
                "output_types": {"rep": "report.html"},
            },
        ],
        "edges": [
            {"from_node": "n1", "from_output": "vcf", "to_node": "n2", "to_input": "input_vcf"},
            {"from_node": "n1", "from_output": "vcf", "to_node": "n3", "to_input": "input_vcf"},
        ],
        "parameter_sweeps": {
            "n2.alpha": [0.1, 0.2],
            "n3.beta": [1, 2],
        },
    }


def test_import_validate_export_and_plan_workflow():
    client = TestClient(app)
    headers = auth_header(client)

    create = client.post("/workflows/import", headers=headers, json=valid_workflow("wf-a"))
    assert create.status_code == 201
    assert create.json()["report"]["acyclic"] is True
    assert "n1" in create.json()["report"]["branching_nodes"]

    listed = client.get("/workflows", headers=headers)
    assert listed.status_code == 200
    assert listed.json()["workflows"][0]["workflow_id"] == "wf-a"

    validate = client.post("/workflows/wf-a/validate", headers=headers)
    assert validate.status_code == 200
    assert validate.json()["report"]["topological_order"][0] == "n1"

    plan = client.get("/workflows/wf-a/plan", headers=headers)
    assert plan.status_code == 200
    assert "n1" in plan.json()["branching_nodes"]

    export = client.get("/workflows/wf-a/export", headers=headers)
    assert export.status_code == 200
    assert export.json()["workflow"]["workflow_id"] == "wf-a"


def test_workflow_rejects_io_type_mismatch():
    client = TestClient(app)
    headers = auth_header(client)
    bad = valid_workflow("wf-bad-io")
    bad["nodes"][2]["input_types"]["input_vcf"] = "reads.fastq.gz"
    response = client.post("/workflows/import", headers=headers, json=bad)
    assert response.status_code == 400
    assert "I/O type mismatch" in response.json()["detail"]


def test_workflow_rejects_cycle():
    client = TestClient(app)
    headers = auth_header(client)
    bad = valid_workflow("wf-cycle")
    bad["edges"].append({"from_node": "n2", "from_output": "rep", "to_node": "n1", "to_input": "vcf"})
    bad["nodes"][0]["input_types"] = {"vcf": "report.html"}
    response = client.post("/workflows/import", headers=headers, json=bad)
    assert response.status_code == 400
    assert "cycle" in response.json()["detail"].lower()


def test_workflow_parameter_sweep_expand():
    client = TestClient(app)
    headers = auth_header(client)
    create = client.post("/workflows/import", headers=headers, json=valid_workflow("wf-sweep"))
    assert create.status_code == 201
    expand = client.post("/workflows/wf-sweep/sweeps/expand", headers=headers)
    assert expand.status_code == 200
    assert expand.json()["count"] == 4
    assert len(expand.json()["expanded_runs"]) == 4


def test_workflow_distributed_execution():
    client = TestClient(app)
    headers = auth_header(client)
    create = client.post("/workflows/import", headers=headers, json=valid_workflow("wf-dist"))
    assert create.status_code == 201
    execute = client.post(
        "/workflows/wf-dist/execute/distributed",
        headers=headers,
        json={"max_workers": 2, "limit_runs": 3},
    )
    assert execute.status_code == 200
    body = execute.json()
    assert body["summary"]["submitted_runs"] == 3
    assert body["summary"]["completed_runs"] == 3
    run_id = body["summary"]["run_id"]

    fetched = client.get(f"/workflows/runs/{run_id}", headers=headers)
    assert fetched.status_code == 200
    assert fetched.json()["summary"]["run_id"] == run_id


def test_workflow_cloud_export_import_local_fallback():
    client = TestClient(app)
    headers = auth_header(client)
    create = client.post("/workflows/import", headers=headers, json=valid_workflow("wf-cloud"))
    assert create.status_code == 201

    exported = client.post("/workflows/wf-cloud/cloud/export", headers=headers)
    assert exported.status_code == 200
    assert exported.json()["workflow_id"] == "wf-cloud"
    from api import workflows

    workflows.reset_workflow_store_for_tests()
    imported = client.post("/workflows/cloud/import/wf-cloud", headers=headers)
    assert imported.status_code == 200
    assert imported.json()["workflow_id"] == "wf-cloud"

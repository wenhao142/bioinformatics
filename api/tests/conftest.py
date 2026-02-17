import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture(autouse=True)
def reset_in_memory_state(monkeypatch):
    # Tests rely on deterministic in-memory state across modules.
    monkeypatch.setenv("VARIANTS_USE_DB", "false")

    from api import audit, datasets, literature, omics, rbac, runs, variants

    if hasattr(variants.STORE, "rows"):
        variants.STORE.rows.clear()
    elif hasattr(variants.STORE, "conn"):
        with variants.STORE.conn.cursor() as cur:
            cur.execute("DELETE FROM variants")

    datasets._DATASETS.clear()
    datasets._NEXT_ID = 1
    omics._EXPR_TABLE.clear()
    literature._PUBMED_RECORDS.clear()
    literature._PUBMED_KEYS.clear()
    rbac.PROJECTS.clear()
    runs._RUNS.clear()
    audit._EVENTS.clear()

    if audit._STORE and hasattr(audit._STORE, "conn"):
        with audit._STORE.conn.cursor() as cur:
            cur.execute("DELETE FROM audit_events")

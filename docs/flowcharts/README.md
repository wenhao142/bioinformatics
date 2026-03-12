# Flowcharts

This folder contains product flowcharts derived from `AGENTS.md`.

Files:

- `LOCAL_ONLY_FLOW.md`
  - Current target flow for local-only execution
  - Local UI + local artifact storage + PostgreSQL metadata
- `FUTURE_CONNECTED_FLOW.md`
  - Future connected flow
  - Local persistence first, then optional cloud dispatch

Rules used to draw these flows:

- Every major step must have validation
- Local persistence happens before cloud dispatch
- Workflow definitions and run metadata must not live only in transient memory
- UI should show where workflows and outputs are stored

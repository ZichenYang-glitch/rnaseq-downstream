# Test suites

The P0 test suite is split by evidence level so that ordinary software tests
cannot be confused with scientific certification.

| Directory | Marker | Purpose | P0 foundation status |
| --- | --- | --- | --- |
| `unit/` | `unit` | Pure Python error and response-contract checks | Active |
| `integration/` | `integration` | Installed CLI behavior across process boundaries | Active |
| `oracle/` | `oracle` | Locked edgeR comparison using the airway fixture | Scaffold only; skipped |
| `simulation/` | `simulation` | Locked compcodeR FDR/TPR gate | Scaffold only; skipped |

Run the active foundation suites after installing the project in editable
mode:

```bash
python -m pip install --no-deps --no-build-isolation -e .
python -m pytest tests/unit tests/integration -v
```

The skipped oracle and simulation modules deliberately fail to confer any
certification. They must only be enabled together with their locked R
implementations, fixtures, tolerances, and machine-readable benchmark
artifacts in later P0 work items.

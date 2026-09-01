# Test suites

The P0 test suite is split by evidence level so that ordinary software tests
cannot be confused with scientific certification.

| Directory | Marker | Purpose | Checkpoint-A status |
| --- | --- | --- | --- |
| `unit/` | `unit` | Contracts, input semantics, provenance, data integrity, and QC math | Active; optional legacy tests may skip when scientific Python packages are absent |
| `integration/` | `integration` | Installed CLI and atomic evidence behavior across process boundaries | Active |
| `oracle/` | `oracle` | Locked edgeR comparison using the airway fixture | Scaffold only; skipped |
| `simulation/` | `simulation` | Locked compcodeR FDR/TPR gate | Scaffold only; skipped |

Run the active standard-library input/CLI suites in the approved P0 environment
after installing the project in editable mode:

```bash
python -m pip install --no-deps --no-build-isolation -e .
python -m pytest tests/unit tests/integration -v
```

The P0 Python lock intentionally has no NumPy, Pandas, PyDESeq2, or plotting
stack. Tests that directly cover the retained experimental data/QC code use
explicit `importorskip` boundaries and must be reported separately with their
Python/scientific-package versions. A larger pass count in that host lane is
expected and is not certification evidence.

The skipped oracle and simulation modules deliberately fail to confer any
certification. They must only be enabled together with their locked R
implementations, fixtures, tolerances, and machine-readable benchmark
artifacts in later P0 work items.

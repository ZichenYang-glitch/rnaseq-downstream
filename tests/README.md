# Test suites

The P0 test suite is split by evidence level so that ordinary software tests
cannot be confused with scientific certification.

| Directory | Marker | Purpose | P0 status |
| --- | --- | --- | --- |
| `unit/` | `unit` | Contracts, input semantics, provenance, result verification, data integrity, and QC math | Active; optional legacy tests may skip when scientific Python packages are absent |
| `integration/` | `integration` | Installed CLI, atomic evidence publication, public input routes, and independent R process behavior | Active; locked-R cases skip only when the runtime is not declared |
| `oracle/` | `oracle` | Same-environment direct edgeR versus toolkit comparison using airway 1.32.0 | Archived evidence always checked; live gate requires locked R |
| `simulation/` | `simulation` | Ten-replicate locked compcodeR 1.48.0 FDR/TPR gate | Archived evidence always checked; live gate requires locked R |

Run the active standard-library input/CLI suites in the approved P0 environment
after installing the project in editable mode:

```bash
python -m pip install --no-deps --no-build-isolation -e .
python -m pytest tests/unit tests/integration -v
```

The locked Python layer includes NumPy 2.4.4 for the C1 display-only PCA path.
It intentionally has no Pandas, PyDESeq2, or plotting stack. Tests that directly
cover retained experimental code with optional dependencies use explicit
`importorskip` boundaries and must be reported separately with their package
versions. A larger pass count in that host lane is expected and is not
certification evidence.

The archived oracle and simulation report tests run even without R. They verify
that each report says `pass`, that every named assertion passed, that the report
contains no machine-local path, and that its implementation inventory still
matches the current statistical code and lock files by SHA-256. Changing those
files therefore invalidates the archived evidence until the locked gates are
rerun.

Ordinary developer runs may omit `RNASEQ_P0_R_LIBRARY`. Only the expensive live
oracle and simulation cases then skip; this is an explicit missing-runtime skip,
not certification evidence. Run the locked lanes with:

```bash
RNASEQ_P0_RSCRIPT=/absolute/path/to/locked/bin/Rscript \
RNASEQ_P0_R_LIBRARY=/absolute/path/to/restored/R/library \
python -m pytest tests/integration tests/oracle tests/simulation -v
```

Certification mode must additionally set
`RNASEQ_P0_REQUIRE_BENCHMARKS=1`. In that mode the oracle/simulation tests fail
if `RNASEQ_P0_R_LIBRARY` is missing, so an absent runtime cannot become a green
skip. The live test locates `Rscript` next to the active locked Python
executable; CI also sets `RNASEQ_P0_RSCRIPT` for the public backend integration
cases. `RNASEQ_P0_BENCHMARK_REPORT_DIR` may point to a non-symlink directory in
which both deterministic live reports are retained for CI archival.

The integration lane also contains a cross-process public-chain test. With
`RNASEQ_P0_RSCRIPT` and `RNASEQ_P0_R_LIBRARY` set it executes
`inspect -> validate -> run -> summarize` on a small eligible featureCounts
fixture, verifies one-JSON stdout at every boundary, and subjects the published
five-file run to the public summary verifier. Certification mode forbids this
test from skipping when either runtime variable is absent.

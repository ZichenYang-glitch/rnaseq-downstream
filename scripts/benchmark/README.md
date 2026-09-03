# P0 certification benchmark runners

These runners create machine-readable evidence for the locked edgeR backend.
They are deliberately separate from the public CLI and must be run with the
same P0 Conda prefix and R library as the toolkit backend.

Both gates in this directory are **backend-kernel** gates. The locked `airway`
object contains integer counts, but this repository does not claim that its
packaged object is accompanied by the toolkit's required featureCounts origin
evidence. The compcodeR matrices are synthetic negative-binomial counts, not
featureCounts outputs. Consequently, neither benchmark fabricates a
featureCounts manifest and neither report certifies the public input route.
They enter the backend through a private benchmark adapter that is unreachable
from the public CLI.

## Airway same-engine oracle

`prepare_airway_fixture.R` extracts the 63,677 by 8 integer assay from the
locked airway 1.32.0 package. `run_airway_direct_oracle.R` independently runs
the literal edgeR 4.10 route:

```text
filterByExpr(design)
normLibSizes(method="TMM")
glmQLFit(legacy=FALSE, robust=TRUE)
glmQLFTest(coef="dextrt")
```

The Python runner compares the direct output with the toolkit backend output.
The tested gene set and coefficient columns must be identical. Every fitted
coefficient, `logFC`, quasi-likelihood `F` statistic, raw `PValue`, and BH
`FDR` must satisfy `math.isclose` with `rel_tol=1e-6` and `abs_tol=1e-10`.
Each result field has an independent assertion and difference metrics in the
machine-readable report. No checked-in expected values are used, so the oracle
cannot pass merely by copying a frozen expected-value file.

## compcodeR negative-binomial gate

The simulation runner generates ten deterministic compcodeR 1.48.0 data sets:
5,000 retained genes, 500 truly differential genes, six samples per condition,
balanced up/down regulation, effect size 2, and no injected outliers. Seeds
5001 through 5010 are fixed. `filter.threshold.total=0` retains all simulated
genes, including all-zero genes, so the truth denominator remains exactly 500;
the toolkit's own `filterByExpr` status remains part of the evaluated route.
The toolkit calls tested genes at BH FDR 0.05.

The CI thresholds are:

- mean replicate false discovery proportion (the Monte Carlo FDR estimate) at
  most 0.065;
- worst-replicate false discovery proportion at most 0.10;
- mean TPR at least 0.50;
- worst-replicate TPR at least 0.45.

These limits were selected after one exploratory locked-runtime precursor run,
which observed mean replicate FDP 0.05289, discovery-weighted pooled FDP
0.05302, worst-replicate FDP 0.08475, mean TPR 0.5717, and worst-replicate TPR
0.524. Review then found that the precursor used
`filter.threshold.total=1`, which removed two all-zero genes (including one
true DE gene) in one replicate. The generator was corrected to retain all
genes, and the already selected limits were left unchanged before the corrected
gate was rerun. They are regression limits, not preregistered evidence from an
untouched evaluation set.

The mean-FDP limit allows finite-simulation variation around the nominal 0.05
target while remaining tight enough to catch a material calibration regression.
The per-replicate FDP guard prevents an acceptable mean from hiding a
catastrophic replicate. The discovery-weighted pooled FDP remains in the report
as a descriptive metric and is not called or gated as FDR. The TPR limits are
regression gates for this frozen scenario, not a universal biological power
claim.

## Reports

After their declared runtime and report destination can be opened, both runners
write a strict JSON report conforming to `benchmark-report-v1.schema.json` for
either a passing gate or an execution/gate failure. Argument-parsing and
unreadable-runtime failures can occur before a report can be initialized. The
report records runtime versions, exact input hashes, thresholds, metrics, and
each gate decision. Its implementation inventory also binds the Conda lock,
renv lock, primary R source lock, top-level environment specification, and
runtime verifier by SHA-256. Reports contain no NaN or Infinity values. A
report with `status: "pass"` is evidence only for its named benchmark and
locked runtime.

The pytest modules under `tests/oracle` and `tests/simulation` invoke these
runners when `RNASEQ_P0_R_LIBRARY` points at the restored locked library. If the
variable is absent, ordinary developer test runs skip the expensive locked
gates. Certification CI must set `RNASEQ_P0_REQUIRE_BENCHMARKS=1` as well as
the library variable; in that mode a missing runtime is a test failure and the
gate cannot silently skip. When `RNASEQ_P0_BENCHMARK_REPORT_DIR` is set, each
live test writes its deterministic report filename into that non-symlink
directory so CI can archive reports even when a gate fails.

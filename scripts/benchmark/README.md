# Locked certification benchmark runners

These runners create machine-readable evidence for the locked edgeR/limma and
DESeq2 statistical backends.
They are deliberately separate from the public CLI and must be run with the
same P0 Conda prefix and R library as the toolkit backend.

The four historical edgeR/limma gates in this directory are
**backend-kernel** gates. The locked
`airway` object contains integer counts, but this repository does not claim that
its packaged object is accompanied by the toolkit's required featureCounts
origin evidence. The compcodeR matrices are synthetic negative-binomial counts,
not featureCounts outputs. Consequently, no benchmark fabricates a featureCounts
manifest and no report certifies the public input route. They enter the backend
through a private benchmark adapter that is unreachable from the public CLI.

The newer DESeq2 runners instead exercise the public CLI with explicitly
self-attested benchmark manifests. This proves numerical routing through
`inspect -> validate -> run -> summarize`, but still does not authenticate the
producer origin of airway serialization or synthetic compcodeR counts.

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

## Airway edgeR/limma pathway same-engine oracle

`run_airway_pathway_oracle.py` reuses the locked airway assay and paired
`~ cell + dex` design. It creates a deterministic local GMT and one-to-one
synthetic symbol annotation with four 50-gene tested sets plus below-minimum and
unmapped controls. Some fixture sets are selected from crude group means to
exercise method output; they are not curated pathways or independent truth.

`run_airway_pathway_direct_oracle.R` is an independent literal implementation
of the same DGEGLM dispatches used by the toolkit: one multi-set mroast call with
9,999 rotations and seed 1729, fry, and camera with per-set inter-gene
correlation estimation. The gate requires exact set/method/hypothesis and status
grids, then compares P values, FDR, mroast proportions, camera raw/effective
correlations, and VIF with `rtol=1e-6` and `atol=1e-10`. It independently
recalculates BH in Python with absolute tolerance `1e-12`.

This oracle detects integration and numerical drift between two independently
written routes through the same locked edgeR/limma engine. It is not external
validation of the methods' statistical assumptions.

## compcodeR pathway-level gate

`generate_pathway_compcoder_fixture.R` and
`run_pathway_compcoder_gate.py` exercise the self-contained directional mroast
and fry families. The runner creates 20 mixed and 40 complete-null fixtures,
each with 5,000 genes and six samples per condition. Mixed fixtures contain 500
true DE genes. Every fixture has 100 null sets (25 each of sizes 20, 40, 80, and
160); mixed fixtures also have five 40-gene Up sets and five 40-gene Down sets.

Set membership is generated from truth categories with a separate RNG stream:
null sets sample truth-null genes without replacement within each set, while
positive sets are disjoint within direction and contain only same-direction
truth-DE genes. Membership never consults observed counts, fitted statistics,
P values, or ranks. Null sets may overlap one another.

At BH FDR 0.05, each of mroast and fry must meet all frozen limits:

- mean mixed-scenario FDP at most 0.10;
- worst mixed-scenario FDP at most 0.25;
- mean mixed-scenario TPR at least 0.80;
- worst mixed-scenario TPR at least 0.60;
- zero wrong-direction significant positive sets;
- no more than four family rejections across 40 complete-null replicates.

The runner also independently verifies BH within every
contrast/method/hypothesis family to absolute tolerance `1e-12`. These limits
were frozen after a disclosed exploratory design pilot. The finite scenarios
are a gross-regression check, not universal evidence of pathway-level FDR
control or biological power. The fixed mroast seed is shared across replicates,
which can induce Monte Carlo dependence and further limits that inference.
Camera is excluded from this self-contained FDR/TPR gate and is covered by the
airway pathway oracle instead.

## Airway DESeq2 Wald/LRT same-engine oracle

`run_deseq2_airway_oracle.py` serializes the locked airway assay as a typed,
self-attested integer-matrix fixture and exercises the complete public CLI for
both Wald and LRT requests. `run_deseq2_airway_direct_oracle.R` independently
runs DESeq2 1.52.0 on the same bytes with the paired `~ cell + dex` design; the
LRT uses reduced `~ cell` and retains one reporting effect without treating it
as the omnibus hypothesis.

The gate requires identical tested-gene sets and compares all fitted
coefficients, log fold changes, statistics, raw P values, and FDR values with
`rtol=1e-6` and `atol=1e-10`. The archived report passes both modes, including
the LRT omnibus P-value/FDR fields. Its scope is same-engine numerical parity
and public routing, not external validation of DESeq2 or proof of featureCounts
origin.

## DESeq2 compcodeR exploration: blocked calibration gate

`run_deseq2_compcoder_gate.py` defines disjoint exploratory seeds 61001--61020
and held-out seeds 62001--62020 for the public DESeq2 Wald route. Before the
exploratory grid was inspected, the runner fixed candidate limits of mean/worst
FDP at most 0.065/0.12 and mean/worst TPR at least 0.45/0.35.

The exploration completed all 20 public CLI chains but observed mean/worst FDP
0.11821/0.13731 at nominal BH FDR 0.05. Direct DESeq2 and independent BH
recalculation reproduced the result. Because this contradicted the expected
conservative behavior and the candidate bounds, the held-out grid was not run,
the limits were not relaxed, and no certification report exists. The archived
exploratory JSON has `evidence_role=threshold_selection_only_not_certification`;
its passing status means that the disclosed execution completed, not that FDR
calibration passed. The detailed audit is in
[`tests/simulation/deseq2-compcoder-method.md`](../../tests/simulation/deseq2-compcoder-method.md).

## Reports

After their declared runtime and report destination can be opened, the runners
write strict JSON reports conforming to
`benchmark-report-v1.schema.json` for either a passing gate or an execution/gate
failure. Argument-parsing and unreadable-runtime failures can occur before a
report can be initialized. The report records runtime versions, exact input
hashes, thresholds, metrics, and each gate decision. Its implementation
inventory also binds the Conda lock, renv lock, primary R source lock, top-level
environment specification, and runtime verifier by SHA-256. Reports contain no
NaN or Infinity values. A report with `status: "pass"` is evidence only for its
named benchmark and locked runtime.

The checked-in report names are:

- `tests/oracle/airway-benchmark-report.json`;
- `tests/oracle/deseq2-airway-benchmark-report.json`;
- `tests/oracle/pathway-airway-benchmark-report.json`;
- `tests/simulation/compcoder-benchmark-report.json`;
- `tests/simulation/deseq2-compcoder-exploratory-report.json`; and
- `tests/simulation/pathway-compcoder-benchmark-report.json`.

The C2 pathway simulation also has a
[human-readable scope report](../../tests/simulation/pathway-compcoder-benchmark-report.md)
covering its independent-gene model and the binomial operating characteristics
of the complete-null cutoff. It is interpretation documentation; the adjacent
JSON remains the authoritative runner-generated evidence.

The pytest modules under `tests/oracle` and `tests/simulation` invoke passing
live gates when `RNASEQ_P0_R_LIBRARY` points at the restored locked library. If
the variable is absent, ordinary developer test runs skip the expensive locked
gates. Certification CI must set `RNASEQ_P0_REQUIRE_BENCHMARKS=1` as well as
the library variable; in that mode a missing runtime is a test failure and a
declared live gate cannot silently skip. The DESeq2 compcodeR module verifies
the archived exploratory failure disclosure but deliberately does not run the
untouched held-out grid. When `RNASEQ_P0_BENCHMARK_REPORT_DIR` is set, each live
test writes its deterministic report filename into that non-symlink directory
so CI can archive reports even when a gate fails. The locked GitHub Actions job
makes every released live gate non-skippable and uploads
`benchmark-results/*.json` as the `p0-benchmark-reports` artifact for 90 days,
including on failed jobs.

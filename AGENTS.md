# Agent Guide: rnaseq-downstream

This file is the entry point for AI agents (Codex, Kimi Code, DeepSeek harnesses,
or any shell-capable agent) operating this toolkit. Read it fully before issuing
any command. It explains how to drive the toolkit, what each statistical method
does, and which claims the evidence does and does not support.

## What this toolkit is

An evidence-oriented bulk RNA-seq downstream toolkit. It deliberately implements
only a small number of analysis paths and refuses unsupported requests. Paths
are listed as evidence-gated only after their applicable gates pass; DESeq2 is
currently implemented but ungated. Every run bundle is bound to hashed inputs
and a locked runtime. Evidence-gated claims additionally bind to passing,
machine-readable benchmark evidence.

**Maturity: research preview.** No output of this toolkit is, by itself, a claim
that a study or conclusion is publication-ready. Benchmark evidence can establish
same-engine numerical fidelity and, only where an applicable simulation gate
exists and passes, calibration in tested scenarios. It does not certify your
study.

## Non-negotiable operating rules

1. **Never bypass a hard fail.** A refusal (`DESIGN_RANK_DEFICIENT`,
   `CONTRAST_NOT_ESTIMABLE`, `INPUT_EVIDENCE_REQUIRED`, ...) is the tool working
   as designed. Do not edit inputs to sneak past a gate; fix the underlying
   problem or report the refusal to the user.
2. **Never modify published output directories.** Evidence bundles are published
   atomically and are never overwritten. To change an analysis, issue a new run.
3. **stdout is one JSON document.** Parse it; do not scrape stderr (logs only).
4. **One backend per run.** Never merge results from two engines into a single
   conclusion ("intersection as truth" is prohibited). The evidence-gated edgeR
   QL path is the primary; the DESeq2 path is a research-preview alternative.
5. **Offline only.** All gene sets, annotation, and references are frozen local
   files. Do not substitute online resources.
6. **Report the evidence tier honestly.** Distinguish `evidence_gated` paths
   from `implemented_ungated` paths whenever you summarize results for a user.

## The five commands

Entry point: `rnaseq-downstream <command>` or `python -m rnaseq_downstream <command>`.

| Command | Purpose |
| --- | --- |
| `capabilities` | Report implemented commands, evidence-gated vs ungated paths, locked runtime identity. Call this first, every session. |
| `inspect` | Read-only audit of declared inputs: provenance, SHA-256, Salmon metadata and inferential-replicate detection. No analysis. |
| `validate` | Validate input semantics and produce an input-validation bundle (the prerequisite for `run`). |
| `run` | Execute an analysis request against a validated input bundle. Publishes an atomic evidence bundle. |
| `summarize` | Independently verify a run bundle (hashes, schemas, status grid, BH recomputation) and summarize outcomes. Never trusts backend-reported values blindly. |

Response envelope (always): `schema_version`, `command`, `status`
(`success`/`error`/`partial`), `data`, `warnings`, `errors`, `artifacts`.
Exit codes: `0` success · `2` invalid request · `3` scientific validation
failure · `4` backend failure · `5` explicit partial run · `70` internal error.

Stable error codes include: `INVALID_REQUEST`, `INPUT_READ_FAILED`,
`INPUT_EVIDENCE_REQUIRED`, `INPUT_INTEGRITY_FAILED`, `COUNT_VALUES_INVALID`,
`SAMPLE_SET_MISMATCH`, `GENE_ID_INVALID`, `ASSAY_PROTOCOL_REQUIRED`,
`SALMON_OFFSET_REQUIRED`, `DESIGN_RANK_DEFICIENT`, `COVARIATE_CONFOUNDED`,
`CONTRAST_NOT_ESTIMABLE`, `BACKEND_FAILED`, `PARTIAL_RUN`,
`FEATURE_NOT_IMPLEMENTED`, `INTERNAL_ERROR`. Treat them as the machine-readable
contract; match on `code`, not on message text.

## Standard workflow

```bash
rnaseq-downstream capabilities
rnaseq-downstream inspect   --request examples/input-requests/salmon-full-length.request.json
rnaseq-downstream validate  --request <input-request.json>
rnaseq-downstream run       --request <analysis-request.json>
rnaseq-downstream summarize <run-bundle-dir>
```

Request templates live in `examples/input-requests/` and
`examples/analysis-requests/`. Exact schemas: `docs/contracts/input-request-v1.md`,
`docs/contracts/analysis-request-v1.md`, `docs/contracts/cli-v1.md`.

## Input semantics (three validated routes)

Input type must be **declared explicitly**; the toolkit never guesses from file
content or dtype.

| Type | Meaning | Key rules |
| --- | --- | --- |
| `featurecounts_integer` | Integer counts from featureCounts | Strict integer grammar; provenance evidence (manifest or per-sample files) required; a bare merged matrix is refused. |
| `salmon_quant_dirs_full_length` | Per-sample Salmon `quant.sf` directories, full-length protocol | tximport with a mandatory transcript-length offset (`countsFromAbundance="no"`). Inferential replicates: 0 → admitted without uncertainty propagation; exactly 1 → refused; ≥2 consistent → propagated as inferential overdispersion by edgeR. |
| `salmon_quant_dirs_three_prime` | Salmon, 3′ tagged protocol | Requires explicit `assay_protocol: three_prime`. Raw counts with `countsFromAbundance="no"` and **without** a length offset (length correction biases 3′ data). Inferential replicates: 0 → admitted without uncertainty propagation; exactly 1 → refused; ≥2 consistent → admitted and recorded but not used for length correction or inferential-overdispersion adjustment. |

For DESeq2, admitted inferential replicates are recorded but not consumed in
either Salmon route. Full-length construction uses
`DESeqDataSetFromTximport`, whose internal conversion uses R `round()`;
three-prime construction explicitly applies R `round()` before
`DESeqDataSetFromMatrix`. Both routes record pre/post matrix hashes, changed-cell
count, maximum absolute change, and per-sample total-count deltas. R uses IEC
60559 ties-to-even rounding; truncation is not permitted.

Internal primary keys are always stable gene IDs; gene symbols are display-only.
Samples must match the metadata exactly — no silent intersection, filling, or
rounding. The DESeq2 Salmon conversions above are explicit and audited.

**Trust boundary:** a typed featureCounts combined-matrix manifest is
self-attested evidence. The gate verifies internal consistency and digest
binding; it cannot cryptographically prove featureCounts produced the matrix.

## Statistical methods, in detail

### 1. Differential expression: edgeR v4 quasi-likelihood (evidence-gated path `edger_ql_p0_v1`)

Pipeline, fixed and recorded in provenance:

```
filterByExpr(design) → normLibSizes(method="TMM")
  → glmQLFit(abundance.trend=TRUE, robust=TRUE,
             winsor.tail.p=c(0.05, 0.1), legacy=FALSE,
             top.proportion=NULL, keep.unit.mat=FALSE)
  → glmQLFTest(contrast=weights)  (lfc_threshold = 0)
  → glmTreat(contrast=weights, lfc=threshold, null="interval")
      (lfc_threshold > 0)
```

- **`filterByExpr`** removes genes too lowly expressed to be testable, using the
  design matrix to account for the experiment's effective replication when
  setting the expression filter. Filtered genes are reported with status
  `filtered`, never silently dropped.
- **TMM normalization** (`normLibSizes`) estimates per-sample composition
  factors so that count comparisons are not confounded by a few highly expressed
  genes or by library size differences.
- **Quasi-likelihood GLM** (`glmQLFit`): RNA-seq counts are overdispersed
  relative to Poisson; the NB dispersion is estimated per gene, and the *QL
  dispersion* models gene-specific residual variability around the NB
  mean-variance relationship and is moderated toward an abundance trend.
  `robust=TRUE` limits the influence of genes with outlying dispersions. This is
  a design rationale for the robust QL route; the current evidence does not
  isolate robustness as the causal explanation for its calibration difference
  from DESeq2 Wald.
- **`glmQLFTest`** is the quasi-likelihood F-test of a contrast — the default
  significance test.
- **`glmTreat(null="interval")`** is the formal effect-threshold test:
  H0 is |log2FC| ≤ threshold. It is *not* a post-hoc |logFC| filter on top of a
  p-value; use it whenever the scientific question is "changed by at least X-fold".
  Note `glmTreat` does not report an F statistic (`statistic_status:
  not_reported`).
- **Design lint (fail-closed, before any fitting):** QR rank of the design
  matrix, positive residual degrees of freedom, detection of complete confounding
  (`limma::nonEstimable` aliases), and contrast estimability (contrast vector in
  the design row space). Violations abort with `DESIGN_RANK_DEFICIENT` or
  `CONTRAST_NOT_ESTIMABLE` and publish nothing.
- **Designs supported:** additive terms only (factors with explicit reference
  level, continuous covariates). Paired/batch designs are expressed as additive
  blocks (e.g., `subject + condition`). Interactions, splines, and random
  effects are out of contract.
- **Contrasts** are arbitrary weight vectors over design-matrix columns, so
  comparisons like `(A+B)/2 − C` are expressible. FDR is Benjamini–Hochberg
  within each contrast and is **recomputed independently** by `summarize`.

Evidence: same-engine airway oracle (zero observed numerical difference for
coefficient/logFC/PValue/F/FDR in the archived run; gate tolerance rtol=1e-6) +
compcodeR negative-binomial simulation gate (FDP/TPR limits). Scope: backend
kernel; input routes are covered by contract plus locked integration tests.

### 2. Gene-set tests: limma fry / mroast / camera (path `edger_ql_p0_v1_limma_gene_sets_v1`)

Run on the same fit, design, and contrast as the DE analysis.

- **`fry`** (primary) and **`mroast`** (corroborative) are **self-contained**
  tests: H0 is "no gene in this set is differentially expressed". They answer
  "does this pathway change at all?" and are valid for paired/blocked designs
  through the design matrix. `fry` is the fast rotation approximation of
  `roast`; `mroast` runs rotation tests over all sets with a recorded seed
  policy.
- **`camera`** is a **competitive** test: H0 is "genes in this set are no more
  DE than the complement", estimating inter-gene correlation per set
  (`inter.gene.cor=NA`, VIF recorded). It answers "is this set enriched relative
  to background?". **camera currently carries same-engine oracle evidence only —
  no simulation FDR/TPR gate — and is labeled `supplementary`.**
- The two test classes have different null hypotheses. Their results, BH
  families, and report sections are kept strictly separate; never present them
  as one ranked list.
- Gene sets come from frozen local GMT files with version + SHA-256 and a
  per-set ID-mapping audit (mapped/unmapped/ambiguous counts, mapping rate,
  tested-universe counts). Sets below the minimum tested-gene count are reported
  as `not_tested` with an explicit reason.
- **Known evidence boundary (disclosed):** the simulation fixture generates
  genes independently, so within-set correlation calibration is *not* covered by
  the gate.

### 3. DESeq2 path (`deseq2_p1_v1_gate_pending`) — research preview, ungated

Implemented with DESeq2 1.52.0 + apeglm: Wald contrasts, LRT omnibus mode
(a proper nested additive reduced design, exactly one reporting contrast,
`lfc_threshold=0`, omnibus-labeled outputs),
formal threshold testing via `results(lfcThreshold=θ, altHypothesis="greaterAbs")`,
and apeglm LFC shrinkage (only a contrast that is exactly one `+1`
non-intercept design coefficient).

**Read this before using it:**

- Same-engine airway oracle: exact numerical parity (Wald and LRT).
- compcodeR calibration gate (6v6, NB simulation): **failed** — mean FDP 0.11821
  vs the predeclared limit 0.065 (native edgeR QL on identical seeds: 0.04834).
- Mechanism diagnostic (archived, hypothesis-generating): false positives
  concentrate in high-dispersion null genes (top true-dispersion quintile
  contributes ~38.9% of FPs); FP genes show depressed final/true dispersion
  ratios. Independent filtering and fit fallback did not explain the result,
  and the gap persisted after aligning tested-gene families. Direct DESeq2
  reproduced the toolkit result, excluding an observed wrapper/routing numeric
  divergence in this diagnostic, not every possible implementation error.
- Interpretation: the wrapper was numerically faithful to official DESeq2 in the
  archived same-engine oracle, but the DESeq2 result was anti-conservative in
  this tested scenario. Use it for exploratory work and literature
  comparability; use the edgeR QL path as the primary for conclusions. Full
  evidence: [exploratory report](tests/simulation/deseq2-compcoder-exploratory-report.json),
  [method audit](tests/simulation/deseq2-compcoder-method.md), and the
  [human-readable](tests/simulation/deseq2-compcoder-mechanism-diagnostic.md) and
  [machine-readable](tests/simulation/deseq2-compcoder-mechanism-diagnostic.json)
  mechanism diagnostics.

### 4. QC mathematics (display layer)

- The public edgeR display PCA consumes the R backend's post-`filterByExpr`,
  post-TMM logCPM matrix. It selects top-variable genes, centers them, and
  **never re-scales** them. The exported logCPM and PCA are display-only and
  never enter inference.
- A separately tested legacy math utility can fit the joint model
  `Y = X_biology·β + X_nuisance·γ + ε` and remove only the nuisance component
  (equivalent to `limma::removeBatchEffect(..., design=X_biology)`). It rejects
  confounded adjustments, but adjusted PCA is not part of the current public
  v1.1/v1.2 display contract.
- Volcano and MA figures render from published result tables; PCA renders from
  the separately exported display-only logCPM matrix. The display layer performs
  no inferential-statistic recomputation (`statistical_recalculation: false` in
  every display manifest).

### 5. Result semantics

Every gene × contrast row carries exactly one status: `filtered`, `not_tested`,
`not_estimable`, `failed`, or `tested`. Numeric fields are empty only in
combination with an explicit status; there are no ambiguous NAs. `summarize`
recomputes BH within each family, verifies every artifact digest and the
five/six-file inventory, and refuses bundles that fail verification.

## Environment and provenance

Locked runtime: two-layer lock (conda-lock for Python/R 4.6.1/toolchain +
renv.lock for Bioconductor 3.23 packages, source archives pinned by SHA-256).
Environment evolution follows the snapshot policy in `environment/README.md`:
old environment files are preserved in append-only snapshots, and every change
re-runs all released live evidence gates with a compatibility report asserting
byte-identical historical numeric artifacts. CI re-runs those live gates from
scratch on each push. The failed DESeq2 calibration study remains archive-only;
CI blocks on its disclosure and integrity checks, not on a rerun or relaxed gate.

## What is NOT available

Interactions (difference-in-differences), random effects / longitudinal models
(dream), voomLmFit sample weights, rank-based exploratory enrichment (fgsea),
TF activity (decoupleR), transcript-level analysis (DTE/DTU), WGCNA,
meta-analysis. Do not improvise these outside the toolkit; request the addition
through the maintainers so it can be certified.

## When something fails

Report the stable error code and the failing stage to the user verbatim. Do not
retry with mutated inputs, do not round counts yourself, do not relax
thresholds. A refusal with a clear code is a successful outcome of the contract.

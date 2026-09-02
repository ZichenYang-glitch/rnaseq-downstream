# Analysis request contract, version 1.0

An analysis request connects one completed, eligible input-validation bundle to
the locked edgeR v4 quasi-likelihood backend. It is strict JSON: all four root
fields are required and additional fields are rejected.

```json
{
  "schema_version": "1.0",
  "validated_input_bundle": "../evidence/sample-set-1",
  "design": {
    "intercept": true,
    "terms": ["batch", "condition"],
    "variables": {
      "batch": {
        "type": "factor",
        "levels": ["batch1", "batch2"]
      },
      "condition": {
        "type": "factor",
        "levels": ["control", "treated"]
      }
    }
  },
  "contrasts": [
    {
      "contrast_id": "treated_vs_control",
      "weights": {"conditiontreated": 1},
      "lfc_threshold": 0
    },
    {
      "contrast_id": "treated_vs_control_min_1",
      "weights": {"conditiontreated": 1},
      "lfc_threshold": 1
    }
  ]
}
```

Relative `validated_input_bundle` paths are resolved from the analysis request,
not from the process working directory. `run` verifies the bundle manifest,
every member digest and size, the shared plan identity, input eligibility, and
the three-prime execution policy before starting R. The backend then consumes
private copies made while hashing the same source streams. An input-only
validation bundle that is ineligible, incomplete, modified, or contains an
unmanifested entry is refused.

## Design

`design` has exactly three fields:

- `intercept` is a boolean.
- `terms` is a non-empty ordered array of unique metadata column names. P0
  accepts additive terms only; interactions, splines, random effects, and raw R
  formulas are outside this contract.
- `variables` describes exactly the names in `terms`. A `continuous` variable
  has only `{"type":"continuous"}`. A `factor` also has an ordered `levels`
  array with at least two unique values. The first factor level is the
  reference level.

Term names use the conservative grammar `[A-Za-z][A-Za-z0-9_.]*`. Metadata
values must be present for every sample. No factor level is guessed or reordered.

R constructs the matrix from this normalized additive design. With an
intercept, its name is `(Intercept)`. A two-level factor named `condition` with
ordered levels `control`, `treated` normally produces the coefficient
`conditiontreated`; a continuous term uses its term name. `run` records the
actual design columns in `analysis.json` and `design.tsv`.

Before fitting, the backend checks the QR rank with the recorded tolerance,
requires positive residual degrees of freedom, and verifies every contrast is
in the design row space. A rank-deficient or completely confounded design
returns `DESIGN_RANK_DEFICIENT`; an unknown, zero, or non-estimable contrast
returns `CONTRAST_NOT_ESTIMABLE`. Neither failure publishes a result directory.

## Contrasts and threshold tests

`contrasts` is a non-empty array. Each item has exactly:

- `contrast_id`: a unique stable identifier;
- `weights`: a non-empty object mapping exact model-matrix column names to
  finite numeric weights; and
- `lfc_threshold`: a finite, non-negative log2-fold-change threshold.

A threshold of zero dispatches to `edgeR::glmQLFTest`. A positive threshold
dispatches to `edgeR::glmTreat(..., null="interval")`. The latter tests the
declared effect-size threshold; it is not a post-hoc `|logFC|` filter.
`glmTreat` does not report an F statistic in this result contract, so its
`statistic` field is empty and the accompanying fields say
`statistic_type: not_reported_by_glmTreat` and
`statistic_status: not_reported`.

## Fixed P0 backend

The only public P0 analysis backend is:

```text
filterByExpr(design)
  -> normLibSizes(method="TMM")
  -> glmQLFit(robust=TRUE, legacy=FALSE)
  -> glmQLFTest or glmTreat(null="interval")
```

Statistical work runs in an independent `Rscript`; Python owns schema
validation, evidence coupling, orchestration, publication, and summarization.
The R process hard-fails unless its R, Bioconductor, edgeR, tximport, and limma
versions exactly match the locked runtime identity.

The three accepted input semantics keep their checkpoint-A meanings:

- `featurecounts_integer` constructs an `edgeR::DGEList` from validated integer
  counts.
- `salmon_quant_dirs_full_length` runs tximport with
  `countsFromAbundance="no"`, retains inferential replicates, and uses
  `edgeR::DGEListFromTximport`; a transcript-length prior offset is mandatory.
  Zero inferential replicates uses `divide=false`; two or more consistent
  replicates per sample uses `divide=true` and requires finite positive edgeR
  overdispersion estimates. Exactly one replicate per sample is rejected before
  analysis.
- `salmon_quant_dirs_three_prime` requires the explicit three-prime declaration
  and uses raw `txi$counts` in `DGEList`, without a gene-length correction or
  offset.

## Published result bundle

A successful `run` atomically publishes exactly five regular files and never
overwrites an existing directory:

- `analysis.json`: normalized scientific provenance, runtime identity, observed
  input route, design lint evidence, contrasts, pipeline, and gene counts;
- `backend_manifest.json`: execution identity plus SHA-256 and byte size for
  every other member;
- `design.tsv`: `sample_id`, `coefficient`, `value` in long form;
- `coefficients.tsv`: `gene_id`, `status`, `coefficient`, `estimate`, `scale` in
  long form; and
- `results.tsv`: one row for every stable `gene_id` and contrast.

`results.tsv` uses this exact header:

```text
gene_id  contrast_id  status  logFC  unshrunk_logFC  logCPM  statistic  statistic_type  statistic_status  PValue  FDR  test_method  lfc_threshold
```

The separators are tabs. `logFC` is log2; fitted coefficients are natural-log
scale. FDR is Benjamini-Hochberg within each contrast. The status vocabulary is
`filtered`, `not_tested`, `not_estimable`, `failed`, and `tested`; empty numeric
fields are interpreted only together with that explicit status. The current
fail-closed backend publishes complete successful runs with `filtered` or
`tested` rows and publishes no bundle after a design, contrast, or backend
failure.

`summarize` independently verifies the exact inventory, manifest/member
digests and sizes, locked identity, strict TSV schemas, complete unique
gene-by-contrast matrix, and method-specific numeric invariants before counting
statuses and FDR-at-0.05 outcomes. It recomputes Benjamini-Hochberg adjustment
within every contrast rather than trusting the reported FDR column. It refuses the private
`backend_kernel_only` artifacts used by certification benchmarks.

## Scope and trust boundary

This is a research-preview evidence-gated path, not a claim that any individual
study or conclusion is publication-ready. In particular, a typed combined
featureCounts manifest is producer-supplied, self-attested evidence. The input
gate can verify its internal consistency and digest binding; it cannot prove
that featureCounts truly produced the matrix. A forged, self-consistent
manifest can therefore pass the input-origin gate. Preserve upstream command,
log, and raw per-sample evidence and perform independent statistical review.

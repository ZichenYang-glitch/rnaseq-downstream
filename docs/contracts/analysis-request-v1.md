# Analysis request contract, versions 1.0 and 1.1

An analysis request connects one completed, eligible input-validation bundle to
the locked edgeR v4 quasi-likelihood backend. Both versions are strict JSON;
additional fields are rejected rather than ignored. Version 1.0 retains its
original four-field schema and has no display request. Version 1.1 requires one
additional `display` object for deterministic post-result SVG output and may
also declare one optional `gene_sets` object for frozen local gene-set tests.

The display request does not change input construction, filtering,
normalization, model fitting, contrasts, P values, or multiple-testing
adjustment. The private R protocol and statistical result schema remain 1.0
when `gene_sets` is absent; a gene-set request selects private schema 1.1 and
adds one manifested inference table. Display artifacts are evidence-bound
presentation outputs, not an additional statistical backend.

## Version 1.0 request

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

Supplying `display` or `gene_sets` with `schema_version: "1.0"` is an error.
Existing version 1.0 requests continue to publish and summarize the original
five-file result directory without a `display/` subdirectory.

## Version 1.1 display request

Version 1.1 keeps `validated_input_bundle`, `design`, and `contrasts` unchanged
and requires this exact additional object:

```json
{
  "schema_version": "1.1",
  "validated_input_bundle": "../evidence/sample-set-1",
  "design": {
    "intercept": true,
    "terms": ["batch", "condition"],
    "variables": {
      "batch": {"type": "factor", "levels": ["batch1", "batch2"]},
      "condition": {"type": "factor", "levels": ["control", "treated"]}
    }
  },
  "contrasts": [
    {
      "contrast_id": "treated_vs_control",
      "weights": {"conditiontreated": 1},
      "lfc_threshold": 0
    }
  ],
  "display": {
    "fdr_threshold": 0.05,
    "pca_top_n": 500,
    "pca_components": [1, 2]
  }
}
```

`display` has exactly three required fields:

- `fdr_threshold` is a finite number in `[0, 1]`. It controls only which tested
  points are highlighted in volcano and MA plots, using the already reported
  per-contrast condition `FDR <= fdr_threshold`. It never filters or changes a
  result row. In particular, it is not combined with a post-hoc `|logFC|`
  cutoff for `glmTreat` results.
- `pca_top_n` is a positive integer. PCA selects at most this many
  positive-variance genes; ties are resolved by stable `gene_id` rather than
  input row order.
- `pca_components` is an ordered array of exactly two distinct, positive,
  one-based component numbers. Their order defines the horizontal and vertical
  axes. A component number cannot exceed `pca_top_n`, and execution also
  rejects components beyond the dimensions supported by the selected matrix
  after centering.

A version 1.1 display request always asks for the fixed C1 set: one volcano and
one MA plot for every declared contrast, plus one sample PCA plot. All plots use
the fixed SVG renderer. Contrast selection, raster formats, themes, arbitrary
axis mappings, configurable gene-ranking methods, and automatic top-gene labels
are not part of version 1.1.

Volcano and MA rendering consume only verified values already present in
`results.tsv`. Filtered rows have no plottable statistics and are not drawn;
their excluded count is recorded. FDR highlighting is a visual category, not a
new differential-expression call.

For PCA, R exports `edgeR::cpm(y, log=TRUE, prior.count=2)` after
`filterByExpr` and TMM normalization. Python selects the requested top
positive-variance genes and applies centered, unscaled PCA: features are
centered but are not divided by their standard deviations. Component signs are
canonicalized deterministically. The PCA is exploratory and display-only. It
does not remove batch effects, enter the edgeR fit, or alter any differential-
expression result. It is a logCPM PCA, not a VST PCA.

The schemas and manifest scope of the five core statistical files remain
unchanged; their bytes can still reflect the analysis-request identity recorded
as provenance. Before the same atomic `run` publication, version 1.1 adds a
`display/` subdirectory with its own manifest, exact source-value evidence, and
SVG members. The display manifest binds the normalized display request, core
source digests, renderer and PCA method identities, per-plot source mappings,
included/excluded row counts, relative member paths, byte sizes, and SHA-256
digests. Publication is deliberately atomic: a requested display failure,
including an invalid PCA component, publishes neither the statistically valid
core nor a partial display directory. `summarize` accepts the original version 1.0
five-file inventory and verifies the additional display inventory when it is
present.

The unchanged core has no field that asserts whether a display was requested.
Consequently, a detached five-file directory cannot distinguish an original
version 1.0 run from a version 1.1 run whose entire `display/` directory was
removed. Preserve the version 1.1 request and successful `run` JSON receipt
when absence detection matters. A present display is still verified against
all five core hashes and reproduced from its bound values. These manifests are
integrity records, not producer signatures, and a newly fabricated,
self-consistent sidecar is outside their trust guarantee.

Every successful `run` JSON receipt records the complete explicit
`glmQLFit` parameter set. Version 1.1 additionally persists that set in the
display manifest. The historical version 1.0 core schema remains unchanged and
therefore carries only its original partial pipeline record; its remaining
defaults are recoverable from the exact locked edgeR identity, not from new
fields retrofitted into old core artifacts.

## Optional version 1.1 gene-set request

Version 1.1 may add this exact object alongside the required `display` object:

```json
{
  "gene_sets": {
    "gmt": {
      "path": "gene-sets/h.all.v2026.1.symbols.gmt",
      "sha256": "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
      "collection": "MSigDB Hallmark",
      "version": "2026.1",
      "identifier_type": "symbol"
    },
    "annotation": {
      "path": "annotation/symbol-to-ensembl.tsv",
      "sha256": "abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789",
      "name": "reviewed symbol-to-Ensembl map",
      "version": "2026-08-01",
      "gene_id_column": "gene_id",
      "symbol_column": "symbol"
    },
    "minimum_tested_genes": 10
  }
}
```

Both source paths must resolve to non-symlink local regular files relative to
the analysis request. URLs and online library names are rejected. Each source
is captured as UTF-8 bytes and must match its declared lowercase SHA-256. At
execution, the same source is captured again while being copied into the
private R workspace; R reads only that private copy. `collection`, `name`, and
both versions are provenance labels supplied by the user, while the digest is
the byte identity the toolkit can verify. The toolkit does not download or
authenticate MSigDB or another provider.

The GMT must contain one unique set ID, one non-empty description, and at least
one symbol per tab-delimited row. Duplicate set IDs, duplicate symbols within a
set, empty fields, and control characters hard-fail the run. The annotation is
a two-column TSV with the exact header `gene_id<TAB>symbol`. Annotation gene
IDs use the same version-stripping policy as the assay and must remain unique
after that policy. A symbol assigned to more than one gene ID is classified as
ambiguous and excluded; it is never expanded silently. Each set records raw and
unique GMT counts, uniquely mapped, ambiguous and unmapped symbol counts,
mapping rate, mapped stable-ID count, post-`filterByExpr` tested count, filtered
count, and count absent from the assay. Mapping rate is uniquely mapped unique
symbols divided by unique GMT symbols.

The backend also records an ordered-membership digest for each test class.
Because the published result bundle binds but does not copy the full GMT and
annotation sources, this digest is execution evidence rather than a
reconstructable mapping proof. `summarize` checks its syntax, eligible set-ID
order, and equality across contrasts; it cannot recompute the member-level
digest from the published bundle alone. Retain the exact digest-bound GMT and
annotation files for a complete mapping audit.

Sets with fewer than `minimum_tested_genes` in the exact post-filter tested
universe remain in `pathway_results.tsv` with `status: not_tested` and a stable
reason. They are not silently dropped or replenished after filtering. A camera
set also remains `not_tested` when it has fewer than two tested genes or covers
the entire tested universe, because a competitive comparison then has no valid
correlation estimate or complement.

All methods consume the same `DGEGLM` returned by the run's single
`glmQLFit`; no logCPM matrix or separately refitted model enters inference:

- `limma::fry` is the primary self-contained count-data test and reports
  directional and mixed hypotheses;
- `limma::mroast` is a corroborative self-contained test with mean set
  statistic, 9,999 rotations, `midp=false`, fixed RNG kind and seed, and
  directional and mixed hypotheses; and
- `limma::camera` is a supplementary competitive test. It estimates
  inter-gene correlation separately for each set with
  `inter.gene.cor=NA`, reports the raw estimate, and records the non-negative
  effective correlation and variance-inflation factor used under
  `allow.neg.cor=false`.

Self-contained and competitive tests answer different questions and are never
combined into one significance list. `fry` and `mroast` test a zero-effect
gene-set null; camera tests enrichment relative to the complement. The
gene-level `lfc_threshold` used by `glmTreat` is not applied to any gene-set
test. For every tested family, BH adjustment is performed separately within
`contrast_id × method_id × hypothesis`. `summarize` independently recomputes
each family from raw P values and rejects pooled or altered FDR values.

See [`examples/analysis-requests/v1.1-pathways.request.json`](../../examples/analysis-requests/v1.1-pathways.request.json)
for a complete schema example with local synthetic assets.

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

## Published core result bundle

A successful `run` always atomically publishes these five core regular files
and never overwrites an existing directory. Version 1.1 additionally publishes
the independently manifested `display/` subdirectory described above. When
`gene_sets` is requested, schema 1.1 also adds `pathway_results.tsv` as a sixth
root inference artifact and requires the display directory to remain present:

- `analysis.json`: normalized scientific provenance, runtime identity, observed
  input route, design lint evidence, contrasts, pipeline, and gene counts;
- `backend_manifest.json`: execution identity plus SHA-256 and byte size for
  every other member;
- `design.tsv`: `sample_id`, `coefficient`, `value` in long form;
- `coefficients.tsv`: `gene_id`, `status`, `coefficient`, `estimate`, `scale` in
  long form; and
- `results.tsv`: one row for every stable `gene_id` and contrast; and
- `pathway_results.tsv`: one row for every declared contrast, gene set,
  method, and applicable directional or mixed hypothesis, including explicit
  mapping, eligibility, method, correlation/rotation, P-value, FDR, and status
  fields.

`results.tsv` uses this exact header:

```text
gene_id  contrast_id  status  logFC  unshrunk_logFC  logCPM  statistic  statistic_type  statistic_status  PValue  FDR  test_method  lfc_threshold
```

When present, `pathway_results.tsv` uses this exact header:

```text
contrast_id  gene_set_id  gene_set_description  method_id  test_class  hypothesis  inference_role  status  status_reason  direction  proportion_down  proportion_up  p_value  fdr  fdr_family_id  gmt_member_count_raw  gmt_symbol_count_unique  mapped_symbol_count_unique  ambiguous_symbol_count_unique  unmapped_symbol_count_unique  mapping_rate  mapped_gene_id_count_unique  tested_gene_count  filtered_gene_count  tested_universe_gene_count  method_ngenes  correlation_status  correlation_estimate_raw  correlation_effective  vif_used  rotation_status
```

The separators in both headers are tabs. `logFC` is log2; fitted coefficients
are natural-log scale. Gene-level FDR is Benjamini-Hochberg within each
contrast. The status vocabulary is
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

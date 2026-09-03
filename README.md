# RNA-seq Downstream Toolkit

> **Status: research preview.** The repository provides one evidence-gated
> bulk RNA-seq differential-expression path. It does not make an individual
> analysis publication-ready; upstream provenance and statistical review are
> still required.

An auditable, non-interactive toolkit that validates declared bulk RNA-seq
inputs, runs a fixed edgeR quasi-likelihood analysis, and publishes
machine-verifiable results through a JSON-only command-line interface.

## Features

- Explicit input semantics: origin and count meaning are never guessed from a
  filename, file extension, or numeric type.
- Stable `gene_id` keys throughout the statistical path; symbols are display
  annotations only.
- Strict handling with no silent sample intersection, missing-value filling,
  count rounding, or gene-symbol aggregation.
- SHA-256-bound input validation and provenance bundles.
- Fail-closed checks for design rank, residual degrees of freedom, complete
  confounding, and contrast estimability.
- A fixed edgeR v4 quasi-likelihood route with TMM normalization.
- Effect-size threshold testing through `glmTreat`, never post-hoc
  `|logFC|` filtering.
- Explicit `filtered` and `tested` statuses in successful result bundles;
  design, contrast, and backend failures publish no partial bundle.
- Stable JSON responses and error codes suitable for shell automation and
  general-purpose coding agents.
- Optional, deterministic SVG volcano, MA, and sample-PCA displays generated
  by the same analysis run and independently checked by `summarize`.
- Locked runtime, same-engine oracle parity, and negative-binomial simulation
  regression gates.

## Supported inputs

| Declared input semantic | Analysis route |
| --- | --- |
| `featurecounts_integer` | Integer counts from original per-sample files or a typed self-attested manifest → `edgeR::DGEList` |
| `salmon_quant_dirs_full_length` | tximport counts plus transcript-length offset → `edgeR::DGEListFromTximport` |
| `salmon_quant_dirs_three_prime` | Raw tximport counts without length correction; requires `assay_protocol: three_prime` |

Full-length Salmon inputs retain inferential replicates when available. Zero
replicates use the ordinary quasi-likelihood route, exactly one replicate is
rejected, and two or more are propagated through the supported uncertainty
route.

Bare merged Salmon gene-count matrices are rejected because they cannot show
whether a transcript-length offset is required. Request templates are under
[`examples/input-requests/`](examples/input-requests/).

## Analysis path

```text
declared input
  → inspect
  → validate and bind provenance
  → design lint
  → filterByExpr
  → TMM normalization
  → glmQLFit
  → glmQLFTest or glmTreat
  → published result bundle
  → summarize and verify
```

P0 supports additive fixed-effect designs composed of declared factor and
continuous variables. Contrasts use explicit weights over model-matrix
columns. A zero effect threshold uses `glmQLFTest`; a positive log2-fold-change
threshold uses `glmTreat(null="interval")`.

## CLI workflow

Install the package from a checkout and restore the
[locked R/Bioconductor runtime](environment/README.md):

```bash
python -m pip install --no-deps --no-build-isolation -e .
```

Then run the non-interactive workflow:

```bash
rnaseq-downstream capabilities

rnaseq-downstream inspect \
  --request input-request.json

rnaseq-downstream validate \
  --request input-request.json \
  --output-dir evidence/sample-set-1

rnaseq-downstream run \
  --request analysis-request.json \
  --output-dir results/sample-set-1 \
  --rscript /path/to/locked/bin/Rscript \
  --r-library /path/to/locked/R/library

rnaseq-downstream summarize \
  --run-dir results/sample-set-1
```

Every command writes one JSON response to stdout, sends diagnostics to stderr,
never prompts, and returns a non-zero exit code with stable structured errors
when it cannot complete safely.

## Results and optional displays

A successful analysis always publishes five core regular files:

| File | Contents |
| --- | --- |
| `analysis.json` | Provenance, runtime identity, design checks, pipeline, and contrasts |
| `backend_manifest.json` | SHA-256 digest and size of each of the other four files |
| `design.tsv` | Model-matrix values by sample and coefficient |
| `coefficients.tsv` | Fitted coefficients and per-gene status |
| `results.tsv` | Per-gene, per-contrast effect, P value, FDR, method, status, and reported statistic or explicit not-reported state |

[Analysis request version 1.1](docs/contracts/analysis-request-v1.md) also adds
a separately manifested `display/` directory in the same atomic publication.
It contains one FDR volcano plot and one MA plot per contrast, plus one sample
PCA plot. The PCA uses edgeR logCPM values after filtering and TMM
normalization, selects up to the requested number of positive-variance genes,
centers them, and does not scale each gene to unit variance. SVG is the only
display format in version 1.1.

The FDR threshold controls point highlighting only. These displays do not
filter result rows, change differential-expression calls, remove batch effects,
or feed values back into the edgeR analysis.

`summarize` verifies the core inventory, hashes, schemas, runtime identity, row
completeness, status vocabulary, and numeric invariants. When `display/` is
present, it additionally verifies the display manifest and reproduces the
source-derived coordinates and SVG content. Existing version 1.0 five-file
directories remain supported. Because the core schemas intentionally contain
no display-presence marker, retain the version 1.1 request and `run` receipt if
you need to detect removal of the entire optional display directory.

## Current scope

The locked airway and compcodeR benchmarks cover the shared edgeR backend
kernel; integration tests separately cover the public input routes. They do
not authenticate producer origin or validate every possible study design. A
combined featureCounts manifest is self-attested: the toolkit verifies and
binds its declared bytes, but cannot prove who produced them.

Interactions, splines, random effects, repeated-measure models, DTU, network
analysis, activity inference, and additional differential-expression backends
are outside P0. The retained PyDESeq2/GSEAPy/Snakemake workflow is
[experimental](legacy/README.md) and is not part of the evidence-gated path.

For schemas, runtime setup, benchmark methods, and development information,
see the [documentation index](docs/README.md).

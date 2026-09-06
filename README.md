# RNA-seq Downstream Toolkit

> **Status: research preview.** The repository provides one evidence-gated
> edgeR differential-expression path and an implemented DESeq2 path whose
> independent certification gates are still pending. It does not make an
> individual analysis publication-ready; upstream provenance and statistical
> review are still required.

An auditable, non-interactive toolkit that validates declared bulk RNA-seq
inputs, runs an explicitly selected locked statistical backend, and publishes
machine-verifiable results through a JSON-only command-line interface.

## Features

- Explicit input semantics: origin and count meaning are never guessed from a
  filename, file extension, or numeric type.
- Stable `gene_id` keys throughout the statistical path; assay symbols are
  display annotations only, while GMT symbols enter pathway tests only through
  the declared frozen symbol-to-`gene_id` map.
- Strict handling with no silent sample intersection, missing-value filling,
  count conversion, or gene-symbol aggregation. DESeq2's required Salmon-count
  rounding is explicit and audited before and after conversion.
- SHA-256-bound input validation and provenance bundles.
- Fail-closed checks for design rank, residual degrees of freedom, complete
  confounding, and contrast estimability.
- A fixed edgeR v4 quasi-likelihood route with TMM normalization.
- Effect-size threshold testing through `glmTreat`, never post-hoc
  `|logFC|` filtering.
- A separate DESeq2 1.52 implementation with Wald tests, formal
  `greaterAbs` LFC-threshold tests, full-versus-reduced LRTs, and optional
  coefficient-only apeglm shrinkage. Its D2 oracle and simulation gates are
  not yet complete.
- Explicit `filtered` and `tested` statuses in successful result bundles;
  design, contrast, and backend failures publish no partial bundle.
- Stable JSON responses and error codes suitable for shell automation and
  general-purpose coding agents.
- Optional, deterministic SVG volcano, MA, and sample-PCA displays generated
  by the same analysis run and independently checked by `summarize`.
- Optional frozen local GMT analysis on the same fitted model: `limma::fry`
  and `limma::mroast` for self-contained tests, plus separately reported
  `limma::camera` competitive tests with mapping audits.
- Locked runtime, same-engine oracle parity, and negative-binomial simulation
  regression gates.

## Supported inputs

| Declared input semantic | edgeR route | DESeq2 route |
| --- | --- | --- |
| `featurecounts_integer` | `edgeR::DGEList` | `DESeqDataSetFromMatrix` |
| `salmon_quant_dirs_full_length` | `DGEListFromTximport`, including the length offset | `DESeqDataSetFromTximport`, including `avgTxLength` normalization |
| `salmon_quant_dirs_three_prime` | Raw tximport counts without length correction | Explicit R `round()` followed by `DESeqDataSetFromMatrix`, without a length offset |

Full-length Salmon inputs retain inferential replicates when available. Zero
replicates use the ordinary route and exactly one replicate is rejected. edgeR
uses two or more replicates for its supported uncertainty route; DESeq2 imports
and verifies them but records that they are not used for inference.

Bare merged Salmon gene-count matrices are rejected because they cannot show
whether a transcript-length offset is required. Request templates are under
[`examples/input-requests/`](examples/input-requests/).
Version 1.1 and 1.2 analysis-request examples are under
[`examples/analysis-requests/`](examples/analysis-requests/).

## Analysis paths

```text
declared input
  → inspect
  → validate and bind provenance
  → design lint
  → explicitly selected backend
      edgeR: filterByExpr → TMM → glmQLFit → glmQLFTest or glmTreat
      DESeq2: DESeq → Wald/greaterAbs or one-contrast full-vs-reduced LRT
              → optional coefficient-only apeglm LFC shrinkage
  → optional edgeR-only fry / mroast / camera gene-set tests
  → published result bundle
  → summarize and verify
```

Both backends accept additive fixed-effect designs composed of declared factor
and continuous variables. Contrasts use explicit weights over model-matrix
columns. DESeq2 LRT requests compare a full design with a proper nested reduced
design, require exactly one reporting contrast, and require an LFC threshold of
zero because their P values are omnibus rather than contrast-specific.

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

## Results, pathways, and optional displays

A successful analysis always publishes five core regular files:

| File | Contents |
| --- | --- |
| `analysis.json` | Provenance, runtime identity, design checks, pipeline, and contrasts |
| `backend_manifest.json` | SHA-256 digest and size of each of the other four files |
| `design.tsv` | Model-matrix values by sample and coefficient |
| `coefficients.tsv` | Fitted coefficients and explicit per-gene status |
| `results.tsv` | Backend-specific, documented per-gene effects, P values, FDR, method, status, reason, and test identity |

[Analysis request versions 1.1 and 1.2](docs/contracts/analysis-request-v1.md)
can add a separately manifested `display/` directory to an edgeR run in the
same atomic publication.
It contains one FDR volcano plot and one MA plot per contrast, plus one sample
PCA plot. The PCA uses edgeR logCPM values after filtering and TMM
normalization, selects up to the requested number of positive-variance genes,
centers them, and does not scale each gene to unit variance. SVG is the only
display format in version 1.1.

When a version 1.1 or 1.2 edgeR request includes frozen local `gene_sets`, the
same atomic run also publishes `pathway_results.tsv`. Every declared set is retained with
its symbol-to-`gene_id` mapping rate and tested-universe size. Sets ineligible
for a declared test are marked `not_tested` with a reason. Self-contained fry/mroast results
and competitive camera results stay separate, and `summarize` independently
checks their distinct BH correction families.

The FDR threshold controls point highlighting only. These displays do not
filter result rows, change differential-expression calls, remove batch effects,
or feed values back into the edgeR analysis.

`summarize` verifies the core inventory, hashes, schemas, runtime identity, row
completeness, status vocabulary, and numeric invariants. For DESeq2 it also
checks independent-filtering BH families, shrinkage separation, rounding
audit structure, and LRT omnibus constraints. When `display/` is
present, it additionally verifies the display manifest and reproduces the
source-derived coordinates and SVG content. Existing version 1.0 five-file
directories remain supported. Because the core schemas intentionally contain
no display-presence marker, retain the version 1.1 request and `run` receipt if
you need to detect removal of the entire optional display directory.

## Current scope

The locked airway and compcodeR benchmarks cover the shared edgeR kernel and
the optional limma gene-set extension; integration tests separately cover the
public input routes. The C2 compcodeR gate evaluates the self-contained fry and
mroast tests; camera currently has same-engine airway parity rather than a
simulation FDR/TPR gate. These checks do not authenticate producer origin or
validate every possible study design. A combined featureCounts manifest is
self-attested: the toolkit verifies and binds its declared bytes, but cannot
prove who produced them.

The DESeq2 D1 implementation is intentionally listed as gate-pending rather
than evidence-gated until its independent airway and compcodeR D2 reports are
complete. Interactions, splines, random effects, repeated-measure models,
GSEA/fgsea, DTU, network analysis, activity inference, and other differential-
expression backends remain outside the evidence-gated path. The retained
PyDESeq2/GSEAPy/Snakemake workflow is
[experimental](legacy/README.md) and is not part of the evidence-gated path.

For schemas, runtime setup, benchmark methods, and development information,
see the [documentation index](docs/README.md).

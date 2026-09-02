# RNA-seq Downstream Toolkit

> **Status: research preview.** One narrow P0 path now connects validated input evidence to the locked edgeR v4 quasi-likelihood backend and machine-verifiable results. Its same-engine oracle and negative-binomial simulation reports certify the shared backend kernel only; locked integration tests separately exercise the public input routes. Neither evidence layer authenticates producer origin or makes an individual analysis or conclusion publication-ready. Independent statistical and upstream-provenance review remains required.

This project is evolving from a Python-led bulk RNA-seq workflow into a small, auditable toolkit with an agent-friendly command-line contract. The immediate goal is one narrow P0 path whose input semantics, statistical execution, failures, and benchmark evidence can all be inspected by software as well as by a person.

The completed P0 vertical slice covers:

- an installable Python package and non-interactive JSON CLI contract;
- stable, machine-readable error handling at the CLI boundary;
- a two-layer environment-locking foundation;
- strict, provenance-backed validation for the three declared P0 input semantics;
- hard design-rank, residual-degrees-of-freedom, and contrast-estimability
  checks;
- the fixed `filterByExpr -> TMM -> glmQLFit -> glmQLFTest/glmTreat` backend in
  an independent R process;
- an exact five-file result bundle plus fail-closed machine summarization;
- same-environment airway oracle and compcodeR FDR/TPR backend-kernel gates;
- corrected stable-ID, integer-count, sample-alignment, PCA, and visualization-only
  covariate-adjustment behavior in the retained experimental workflow;
- unit, integration, oracle, and simulation test entry points.

This is deliberately one evidence-gated research-preview path, not a
publication-grade end-to-end claim or a general claim of scientific correctness
for arbitrary study designs, upstream data, or interpretation.

## Project Boundaries

The repository currently contains two clearly separate tracks.

### P0 evidence-gated path

The implemented path is:

```text
declared input semantics -> tximport where appropriate -> edgeR v4 QL -> structured results and evidence
```

The input gate accepts only explicitly declared, provenance-backed inputs and
uses stable gene IDs internally. `run` verifies the completed input bundle and
its digests, executes its fixed semantic route, requires a full-rank design
with positive residual degrees of freedom, and refuses non-estimable
contrasts. A zero effect threshold uses `glmQLFTest`; a positive log2-fold
threshold uses `glmTreat(null="interval")`, never post-hoc `|logFC|` filtering.

`summarize` does not trust file presence. It independently verifies the exact
result inventory, manifest digests and sizes, locked backend/runtime identity,
strict TSV schemas, complete gene-by-contrast rows, explicit status vocabulary,
and method-specific numeric invariants before reporting compact counts.

A successful input-only `validate` still does not prove that a design is legal
or that analysis ran. Conversely, a completed run cannot prove producer origin:
the typed combined featureCounts manifest is self-attested. The validator can
bind and check the declared bytes but cannot prove that featureCounts truly
created them. Preserve original per-sample outputs, commands, and logs.

### Retained experimental legacy workflow

The existing `PyDESeq2`/`GSEAPy`/Snakemake workflow remains in the repository so its behavior can be reviewed and migrated incrementally. It is experimental and is outside the P0 evidence-gated DAG.

Checkpoint A removed the known silent missing-value fill, rounding, display-symbol aggregation, sample-set intersection, scaled PCA, and incorrect nuisance-only residualization from this retained path. PCA now uses the top positive-variance genes after the selected transform, centers without feature scaling, and records method/feature evidence. Adjusted PCA uses a joint biology-plus-nuisance model and refuses completely confounded designs. These are bug fixes to an **experimental** PyDESeq2 workflow, not certification of its differential-expression results. The historical HOMER helper was removed from the core CLI and Snakemake DAG and retained only under `extras/homer/` as unsupported reference code.

Corrected PCA and adjusted-PCA values intentionally differ from legacy output. The workflow refuses to overwrite legacy PCA artifacts that lack current method evidence; archive old result directories before rerunning. The duplicated `run_integrated_pydeseq2.py` entry point is fail-closed because it bypassed these corrections.

The legacy config still names `salmon.merged.gene_counts.tsv` as a compatibility placeholder. It now hard-fails if that file contains fractional estimated counts, and the P0 input gate rejects it regardless because a bare merged matrix cannot prove whether a transcript-length offset is required.

## Install the Toolkit Package

For development from a repository checkout:

```bash
python -m pip install -e .
```

The console command is then available as either:

```bash
rnaseq-downstream capabilities
python -m rnaseq_downstream capabilities
```

The package includes the independent R backend script. The exact locked R and
Bioconductor library is still required at execution time; package installation
does not make the legacy workflow part of the evidence-gated path.

## CLI Contract

The command surface is:

```text
rnaseq-downstream capabilities
rnaseq-downstream inspect --request REQUEST.json
rnaseq-downstream validate --request REQUEST.json --output-dir EVIDENCE_DIR
rnaseq-downstream run --request ANALYSIS_REQUEST.json --output-dir RUN_DIR
rnaseq-downstream summarize --run-dir RUN_DIR
```

The contract is intentionally suitable for shell scripts and general-purpose coding agents:

- stdout contains exactly one JSON envelope and no human-oriented logging;
- diagnostics and logs go to stderr;
- commands never prompt interactively;
- every response includes `schema_version`, `command`, `status`, `data`, `warnings`, `errors`, and `artifacts`;
- failures use stable error codes and non-zero process exit codes; and
- help and invalid-argument responses are JSON rather than argparse prose.

The normative field and exit-code definitions are in the [CLI JSON contract](docs/contracts/cli-v1.md).

All five commands are operational. `inspect` is read-only. `validate` performs
strict numeric checks and publishes `validated_request.json`, `input_plan.json`,
`provenance.json`, and a digest-bearing `bundle_manifest.json` as one
non-overwriting bundle. These two input commands explicitly report
`design: "not_run"`, `backend: "not_run"`, `runnable: false`, and
`analysis_path_certified: false`; only `run` evaluates the analysis design and
backend.

The three accepted and executed semantics are:

- `featurecounts_integer`, backed by original per-sample files or a typed,
  digest-bound combined-matrix manifest;
- `salmon_quant_dirs_full_length`, using tximport original counts plus the
  required length offset and `edgeR::DGEListFromTximport`; and
- `salmon_quant_dirs_three_prime`, requiring `assay_protocol: three_prime` and
  using raw `txi$counts` without a length offset.

The normative schemas and examples are in the [input request
contract](docs/contracts/input-request-v1.md) and [analysis request
contract](docs/contracts/analysis-request-v1.md). A bare merged Salmon
gene-count matrix is rejected; input origin is never guessed from a filename or
dtype.

The complete non-interactive flow is:

```bash
rnaseq-downstream inspect --request input-request.json
rnaseq-downstream validate \
  --request input-request.json \
  --output-dir evidence/sample-set-1
rnaseq-downstream run \
  --request analysis-request.json \
  --output-dir results/sample-set-1 \
  --rscript /path/to/locked/bin/Rscript \
  --r-library /path/to/locked/R/library
rnaseq-downstream summarize --run-dir results/sample-set-1
```

The analysis request references the validation bundle and declares factor or
continuous additive terms plus explicit contrast weights keyed by exact design
columns. See the contract before constructing requests programmatically.
Invalid CLI arguments exit `2`; scientific design or contrast failures exit
`3`; backend failures exit `4`. No failed design or backend run publishes a
result bundle.

## Reproducible Environment Strategy

P0 uses an approved two-layer lock because the requested Bioconductor 3.23 packages are not yet available as a complete Bioconda solve:

1. Conda defines and locks the Linux runtime foundation, including Python, R 4.6.x, native libraries, and build/runtime tooling.
2. `renv` locks the R and Bioconductor package graph, including Bioconductor 3.23, edgeR 4.10.0, tximport 1.40.0, limma, airway, and compcodeR.

Any archived R package sources used by the lock are accompanied by a SHA-256 manifest so source identity can be checked independently. The Conda lock, `renv.lock`, and checksum manifest are complementary parts of one environment contract; using only one layer is not the P0 environment.

The legacy [`environment.yaml`](environment.yaml) remains for the experimental Snakemake workflow. It is not the P0 lock and should not be cited as reproducibility evidence for the P0 edgeR path.

## Test Layout

The test suite is divided by evidence type:

```text
tests/
├── unit/         # pure Python contracts and error mapping
├── integration/  # installed CLI and process-boundary behavior
├── oracle/       # same-environment airway edgeR parity gate
└── simulation/   # locked compcodeR FDR/TPR gate
```

Unit and integration tests exercise the contracts, input semantics, public
backend routes, result verifier, and process boundaries. Oracle and simulation
tests also require the exact locked R runtime. Their archived reports are
machine-readable evidence; CI must rerun the live gates rather than treating an
ordinary missing-runtime skip as a pass.

Run the software checks from the locked P0 Conda environment with:

```bash
python -m pip install --no-deps --no-build-isolation -e .
python -m pytest tests/unit tests/integration
```

Run the live certification lanes with the explicit locked runtime:

```bash
RNASEQ_P0_RSCRIPT=/path/to/locked/bin/Rscript \
RNASEQ_P0_R_LIBRARY=/path/to/locked/R/library \
python -m pytest tests/oracle tests/simulation -v
```

The approved P0 lock is intentionally toolchain-only on the Python side and does
not add NumPy/Pandas/PyDESeq2 for the retained legacy workflow. NumPy-only QC
math tests therefore run in a separately reported scientific-dependency lane;
they are not evidence that the legacy path is certified or exactly locked.

## Legacy Workflow Reference

The following interface describes the retained experimental workflow only.

### Components

- [`Snakefile`](Snakefile) and [`workflow/rules/rnaseq.smk`](workflow/rules/rnaseq.smk) orchestrate legacy steps.
- [`workflow_config.yaml`](workflow_config.yaml) configures merged matrices, metadata, contrasts, QC, enrichment, and reports.
- [`modules/`](modules) contains the legacy Python analysis implementation.
- [`scripts/`](scripts) contains the legacy Snakemake entry points.

### Legacy setup

Create or reuse a Snakemake environment, then allow Snakemake to build the rule environment from the legacy `environment.yaml`:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda run -n snakemake snakemake --use-conda --cores 4
```

A dry run is available with:

```bash
conda run -n snakemake snakemake -n --use-conda --cores 4
```

The legacy direct entry point is also retained:

```bash
python main.py
python main.py --step qc
python main.py --step deseq
python main.py --step gsea
python main.py --step report
```

These commands are provided for compatibility and migration work. They do not use the P0 edgeR backend or satisfy its evidence gates.
This checkpoint retains that execution surface as-is and does not validate it.

## Evidence maintenance and future scope

The first narrow path is present, but the repository remains a research
preview. Releasing a change to the statistical kernel requires a fresh locked
airway parity report and compcodeR FDR/TPR report. New backends, repeated-measure
models, DTU, network analysis, and activity inference remain outside P0; they
must earn their own contracts, tests, and benchmark gates rather than inherit
this path's evidence.

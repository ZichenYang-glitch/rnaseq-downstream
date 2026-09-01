# RNA-seq Downstream Toolkit

> **Status: research preview.** Input semantics can now be inspected and validated, but no statistical analysis path has completed the planned design, backend, oracle, or simulation gates. Do not use the current outputs as a substitute for independent statistical review.

This project is evolving from a Python-led bulk RNA-seq workflow into a small, auditable toolkit with an agent-friendly command-line contract. The immediate goal is one narrow P0 path whose input semantics, statistical execution, failures, and benchmark evidence can all be inspected by software as well as by a person.

The completed foundation and checkpoint-A work covers:

- an installable Python package and non-interactive JSON CLI contract;
- stable, machine-readable error handling at the CLI boundary;
- a two-layer environment-locking foundation; and
- strict, provenance-backed validation for the three declared P0 input semantics;
- corrected stable-ID, integer-count, sample-alignment, PCA, and visualization-only
  covariate-adjustment behavior in the retained experimental workflow; and
- unit, integration, oracle, and simulation test entry points.

It does **not** implement or certify the planned edgeR analysis path yet.

## Project Boundaries

The repository currently contains two clearly separate tracks.

### Planned P0 evidence-gated path

The planned path is:

```text
declared input semantics -> tximport where appropriate -> edgeR v4 QL -> structured results and evidence
```

The input gate now accepts only explicitly declared, provenance-backed inputs and uses stable gene IDs internally. It records the required future route, including transcript-length-offset semantics, but does not yet execute that route. Design and contrast estimability checks, edgeR QL fitting, threshold tests, per-gene test states, and benchmark gates belong to checkpoint B. None of those scientific capabilities should be inferred from a successful input-only validation.

### Retained experimental legacy workflow

The existing `PyDESeq2`/`GSEAPy`/Snakemake workflow remains in the repository so its behavior can be reviewed and migrated incrementally. It is experimental and is outside the P0 evidence-gated DAG.

Checkpoint A removed the known silent missing-value fill, rounding, display-symbol aggregation, sample-set intersection, scaled PCA, and incorrect nuisance-only residualization from this retained path. PCA now uses the top positive-variance genes after the selected transform, centers without feature scaling, and records method/feature evidence. Adjusted PCA uses a joint biology-plus-nuisance model and refuses completely confounded designs. These are bug fixes to an **experimental** PyDESeq2 workflow, not certification of its differential-expression results. The historical HOMER helper was removed from the core CLI and Snakemake DAG and retained only under `extras/homer/` as unsupported reference code.

Corrected PCA and adjusted-PCA values intentionally differ from legacy output. The workflow refuses to overwrite legacy PCA artifacts that lack current method evidence; archive old result directories before rerunning. The duplicated `run_integrated_pydeseq2.py` entry point is fail-closed because it bypassed these corrections.

The legacy config still names `salmon.merged.gene_counts.tsv` as a compatibility placeholder. It now hard-fails if that file contains fractional estimated counts, and the P0 input gate rejects it regardless because a bare merged matrix cannot prove whether a transcript-length offset is required.

## Install the Foundation Package

For development from a repository checkout:

```bash
python -m pip install -e .
```

The console command is then available as either:

```bash
rnaseq-downstream capabilities
python -m rnaseq_downstream capabilities
```

The package installation exposes the machine contract; it does not make the legacy scientific workflow part of the planned evidence-gated path.

## CLI Contract

The command surface is:

```text
rnaseq-downstream capabilities
rnaseq-downstream inspect --request REQUEST.json
rnaseq-downstream validate --request REQUEST.json --output-dir EVIDENCE_DIR
rnaseq-downstream run
rnaseq-downstream summarize
```

The contract is intentionally suitable for shell scripts and general-purpose coding agents:

- stdout contains exactly one JSON envelope and no human-oriented logging;
- diagnostics and logs go to stderr;
- commands never prompt interactively;
- every response includes `schema_version`, `command`, `status`, `data`, `warnings`, `errors`, and `artifacts`;
- failures use stable error codes and non-zero process exit codes; and
- help and invalid-argument responses are JSON rather than argparse prose.

The normative field and exit-code definitions are in the [CLI JSON contract](docs/contracts/cli-v1.md).

`capabilities`, `inspect`, and input-only `validate` are operational. `inspect` is read-only. `validate` performs the strict numeric checks and publishes `validated_request.json`, `input_plan.json`, `provenance.json`, and a digest-bearing `bundle_manifest.json` as one non-overwriting bundle. Both commands explicitly report `design: "not_run"`, `backend: "not_run"`, `runnable: false`, and `analysis_path_certified: false`.

The three accepted semantics are:

- `featurecounts_integer`, backed by original per-sample files or a typed,
  digest-bound combined-matrix manifest;
- `salmon_quant_dirs_full_length`, planned for tximport original counts plus
  length offset and `edgeR::DGEListFromTximport`; and
- `salmon_quant_dirs_three_prime`, requiring `assay_protocol: three_prime` and
  planned to use raw `txi$counts` without a length offset.

The normative request schema and templates are in the [input request contract](docs/contracts/input-request-v1.md). A bare merged Salmon gene-count matrix is rejected; input origin is never guessed from a filename or dtype.

`run` and `summarize` remain reserved and deliberately fail with:

```json
{
  "status": "error",
  "errors": [
    {
      "code": "FEATURE_NOT_IMPLEMENTED"
    }
  ]
}
```

The actual response contains the complete envelope described above. These reserved commands exit `2`; they must never be treated as successful analysis. Invalid CLI arguments likewise exit `2`, using `INVALID_REQUEST`. Scientific input-validation failures exit `3` with a specific structured error code.

## Reproducible Environment Strategy

P0 uses an approved two-layer lock because the requested Bioconductor 3.23 packages are not yet available as a complete Bioconda solve:

1. Conda defines and locks the Linux runtime foundation, including Python, R 4.6.x, native libraries, and build/runtime tooling.
2. `renv` locks the R and Bioconductor package graph, including Bioconductor 3.23, edgeR 4.10.0, tximport 1.40.0, limma, airway, and compcodeR.

Any archived R package sources used by the lock are accompanied by a SHA-256 manifest so source identity can be checked independently. The Conda lock, `renv.lock`, and checksum manifest are complementary parts of one environment contract; using only one layer is not the P0 environment.

The legacy [`environment.yaml`](environment.yaml) remains for the experimental Snakemake workflow. It is not the P0 lock and should not be cited as reproducibility evidence for the planned edgeR path.

## Test Layout

The test suite is divided by evidence type:

```text
tests/
├── unit/         # pure Python contracts and error mapping
├── integration/  # installed CLI and process-boundary behavior
├── oracle/       # same-environment edgeR parity fixtures (future P0 work)
└── simulation/   # FDR/TPR gates (future P0 work)
```

At this checkpoint the unit and integration layers exercise the foundation, input semantics, and process boundary. Oracle and simulation entry points are scaffolds only; their existence does not mean the scientific gates have run or passed.

Run the available foundation checks from the locked P0 Conda environment with:

```bash
python -m pip install --no-deps --no-build-isolation -e .
python -m pytest tests/unit tests/integration
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

These commands are provided for compatibility and migration work. They do not use the planned edgeR backend or satisfy the P0 completion criteria.
This checkpoint retains that execution surface as-is and does not validate it.

## Roadmap to the First Evidence-Gated Path

Checkpoint B will add, in order:

1. an independent Rscript edgeR v4 QL backend with hard design-lint failures;
   and
2. same-environment oracle parity plus negative-binomial FDR/TPR simulation
   gates and machine-readable benchmark evidence.

Until those gates are implemented, executed, and archived as machine-readable benchmark reports, the repository remains a research preview.

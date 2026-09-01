# RNA-seq Downstream Toolkit

> **Status: research preview.** The repository is currently building its engineering and reproducibility foundation. No statistical analysis path in this checkpoint has completed the planned evidence gates. Do not use the current outputs as a substitute for independent statistical review.

This project is evolving from a Python-led bulk RNA-seq workflow into a small, auditable toolkit with an agent-friendly command-line contract. The immediate goal is one narrow P0 path whose input semantics, statistical execution, failures, and benchmark evidence can all be inspected by software as well as by a person.

This checkpoint covers only:

- an installable Python package and non-interactive JSON CLI contract;
- stable, machine-readable error handling at the CLI boundary;
- a two-layer environment-locking foundation; and
- unit, integration, oracle, and simulation test entry points.

It does **not** implement or certify the planned edgeR analysis path yet.

## Project Boundaries

The repository currently contains two clearly separate tracks.

### Planned P0 evidence-gated path

The planned path is:

```text
declared input semantics -> tximport where appropriate -> edgeR v4 QL -> structured results and evidence
```

It will accept only explicitly declared, provenance-backed inputs and will use stable gene IDs internally. Design and contrast estimability checks, transcript-length offsets, edgeR QL fitting, threshold tests, per-gene test states, and benchmark gates belong to later P0 work items. None of those scientific capabilities should be inferred from the presence of the current CLI commands or environment files.

### Retained experimental legacy workflow

The existing `PyDESeq2`/`GSEAPy`/Snakemake workflow remains in the repository so its behavior can be reviewed and migrated incrementally. It is experimental and is outside the P0 evidence-gated DAG.

In particular, the legacy path still starts from merged gene-level matrices and includes known input-semantics, identifier, PCA, and covariate-adjustment issues scheduled for later work. Its reports and statistical outputs are not evidence that the planned P0 path has passed validation. Optional HOMER support is also legacy functionality, not a core P0 capability.

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

The reserved command surface is:

```text
rnaseq-downstream capabilities
rnaseq-downstream inspect
rnaseq-downstream validate
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

At this checkpoint, only `capabilities` is operational. It exits `0` with `status: "success"` and describes the current surface. `inspect`, `validate`, `run`, and `summarize` are reserved for subsequent P0 work and deliberately fail with:

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

The actual response contains the complete envelope described above. These reserved commands exit `2`; they must never be treated as successful analysis. Invalid CLI arguments likewise exit `2`, using `INVALID_REQUEST`.

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

At this checkpoint the unit and integration layers exercise the foundation. Oracle and simulation entry points are scaffolds only; their existence does not mean the scientific gates have run or passed.

Run the available foundation checks from the locked P0 Conda environment with:

```bash
python -m pip install --no-deps --no-build-isolation -e .
python -m pytest tests/unit tests/integration
```

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
python main.py --step motif
python main.py --step report
```

These commands are provided for compatibility and migration work. They do not use the planned edgeR backend or satisfy the P0 completion criteria.
This checkpoint retains that execution surface as-is and does not validate it.

## Roadmap to the First Evidence-Gated Path

Subsequent P0 work will add, in order:

1. corrections for known legacy data handling and QC errors;
2. explicit, provenance-backed input semantics;
3. an independent Rscript edgeR QL backend with hard design-lint failures; and
4. same-environment oracle parity plus negative-binomial FDR/TPR simulation gates.

Until those gates are implemented, executed, and archived as machine-readable benchmark reports, the repository remains a research preview.

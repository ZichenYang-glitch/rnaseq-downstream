# Experimental legacy workflow

This directory contains the retained PyDESeq2, GSEAPy, and Snakemake workflow.
It is available for compatibility and migration work, but it is not part of
the P0 evidence-gated edgeR path.

## Contents

- `main.py` provides the monolithic experimental runner.
- `modules/` contains its data, PyDESeq2, enrichment, and report code.
- `scripts/` contains step-specific Snakemake entry points.
- `workflow/rules/rnaseq.smk` defines the experimental DAG.
- `workflow_config.yaml` and `environment.yaml` configure that DAG.
- `extras/homer/` preserves the unsupported historical HOMER helper.
- `examples/` contains optional legacy manifest and contrast templates.

## Run

From the repository root, create a Snakemake environment and let Snakemake
create the rule environment:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda run -n snakemake snakemake --use-conda --cores 4
```

The root `Snakefile` delegates to `legacy/Snakefile`, so the ordinary Snakemake
command remains available. To run the monolithic interface instead:

```bash
RNASEQ_CONFIG=legacy/workflow_config.yaml python -m legacy
RNASEQ_CONFIG=legacy/workflow_config.yaml python -m legacy --step qc
```

The step modules use the same checkout-root module form, for example:

```bash
RNASEQ_CONFIG=legacy/workflow_config.yaml python -m legacy.scripts.run_qc
```

## Boundaries

The default legacy configuration still names a merged Salmon gene-count matrix
for historical compatibility. Fractional estimated counts hard-fail, and the
P0 input gate rejects bare merged matrices regardless. Do not interpret this
legacy route as equivalent to the P0 tximport and edgeR path.

PCA uses top positive-variance genes after transformation and centers without
feature scaling. Adjusted PCA removes nuisance terms while retaining declared
biology and refuses completely confounded designs. These are integrity fixes
for visualization outputs; they do not certify the legacy differential-
expression backend.

Corrected PCA values intentionally differ from older runs. Archive existing
result directories before rerunning. The duplicated
`run_integrated_pydeseq2.py` entry point remains fail-closed because it bypassed
the corrected implementation.

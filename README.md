# RNA-seq Downstream Workflow

This repository contains a modular RNA-seq downstream workflow built around `PyDESeq2`, `GSEAPy`, optional `HOMER`, and a Snakemake execution layer.

The current recommended mode is:

```bash
conda run -n snakemake snakemake --use-conda --cores 4
```

In this model:

- your existing `snakemake` conda environment provides the `snakemake` executable
- this repository's [environment.yaml](/home/irenadler/projects/rnaseq-downstream/environment.yaml) provides the rule runtime dependencies
- if you do not pass `--conda-prefix`, Snakemake will cache rule environments under `.snakemake/conda/` inside this repository

## What It Does

The workflow currently supports:

- loading merged gene-level count and TPM matrices
- validating metadata, design settings, and contrasts before model fitting
- QC with variance-stabilized counts by default
- differential expression with `PyDESeq2`
- per-contrast volcano plots
- preranked GSEA
- optional motif analysis through `HOMER`
- an integrated expression/statistics summary table

## Workflow Layout

```text
rnaseq-downstream/
├── Snakefile
├── workflow_config.yaml
├── environment.yaml
├── main.py
├── scripts/
│   ├── run_qc.py
│   ├── run_deseq.py
│   ├── run_gsea.py
│   ├── run_report.py
│   └── run_motif.py
├── workflow/
│   └── rules/
│       └── rnaseq.smk
└── modules/
    ├── data.py
    ├── deseq.py
    ├── enrichment.py
    ├── motif.py
    └── report.py
```

## Inputs

You need three main input files:

1. count matrix
   gene-level raw counts, tab-separated
2. TPM matrix
   gene-level TPM matrix, tab-separated
3. metadata table
   one row per sample, including a sample ID column and the design column

Expected metadata format:

```tsv
sample_id	group
sample_1	M0
sample_2	M1
sample_3	M2
```

Count and TPM matrices should contain genes in rows and sample IDs in columns.

## Configuration

Main runtime settings live in [workflow_config.yaml](/home/irenadler/projects/rnaseq-downstream/workflow_config.yaml).

Important fields:

- `COUNTS_FILE`: merged gene count matrix
- `TPM_FILE`: merged gene TPM matrix
- `METADATA_FILE`: sample metadata
- `OUTPUT_DIR`: output root directory
- `ANNOTATION_FILE`: optional local annotation table for gene ID to gene name mapping
- `ANNOTATION_GENE_ID_COL`: gene ID column in the annotation table
- `ANNOTATION_GENE_NAME_COL`: gene name column in the annotation table
- `DESIGN_FACTOR`: metadata column used in DESeq2
- `DESIGN`: formula passed to `PyDESeq2`, for example `~ group` or `~ batch + group`
- `REFERENCE_LEVEL`: baseline group
- `REFERENCE_LEVELS`: mapping used to force stable categorical baselines
- `CONTINUOUS_FACTORS`: metadata columns to treat as numeric covariates
- `CONTRASTS_FILE`: optional contrast table path
- `CONTRASTS`: list of `[treatment, control]`
- `MIN_COUNTS`: low-count gene filter threshold
- `PADJ_THRESH`: adjusted p-value cutoff
- `LOGFC_THRESH`: log2 fold-change cutoff
- `QC_TRANSFORM`: `vst` or `log1p`
- `QC_ADJUST_FACTORS`: nuisance covariates to regress out before adjusted PCA
- `QC_ANNOTATION_FACTORS`: metadata fields used as clustermap annotation color bars
- `VST_USE_DESIGN`: whether VST should use the design matrix
- `GSEA_GENE_SETS`: gene set libraries passed to `GSEAPy`
- `GSEA_PERMUTATIONS`: permutation count for preranked GSEA
- `GSEA_RANK_METRIC`: defaults to `stat`, falls back to `log2FoldChange` if absent
- `SHRINK_LFC`: whether to apply LFC shrinkage where safe
- `SHRINK_LFC_ADAPT`: whether shrinkage prior should adapt to the data
- `RUN_MOTIF`: enable or disable HOMER
- `HOMER_SPECIES`: species passed to HOMER
- `N_CPUS`: CPUs used inside PyDESeq2 and related steps

Example:

```yaml
COUNTS_FILE: salmon.merged.gene_counts.tsv
TPM_FILE: salmon.merged.gene_tpm.tsv
METADATA_FILE: metadata.txt
OUTPUT_DIR: Final_Analysis_Pipeline

DESIGN_FACTOR: group
REFERENCE_LEVEL: M0
REFERENCE_LEVELS:
  group: M0
CONTINUOUS_FACTORS: []
CONTRASTS_FILE: null

CONTRASTS:
  - [M1, M0]
  - [M2, M0]
  - [M1, M2]

QC_TRANSFORM: vst
RUN_GSEA: true
RUN_MOTIF: false
N_CPUS: 4
```

If you prefer a contrast table, point `CONTRASTS_FILE` to a TSV like [contrasts.example.tsv](/home/irenadler/projects/rnaseq-downstream/contrasts.example.tsv):

```tsv
factor	treatment	control	name
group	M1	M0	M1_vs_M0
group	M2	M0	M2_vs_M0
group	M1	M2	M1_vs_M2
```

For a multi-factor design, a common pattern is:

```yaml
DESIGN: ~ batch + group
DESIGN_FACTOR: group
REFERENCE_LEVELS:
  group: M0
  batch: batch1
CONTINUOUS_FACTORS: []
```

## Running The Workflow

### Recommended: Snakemake

Run the full workflow:

```bash
conda run -n snakemake snakemake --use-conda --cores 4
```

This uses the default Snakemake conda cache location:

```text
.snakemake/conda/
```

Use a custom prefix only if you explicitly want a shared cache across projects.

Dry-run:

```bash
conda run -n snakemake snakemake -n --use-conda --cores 4
```

Use a different config file:

```bash
conda run -n snakemake snakemake --use-conda --cores 4 --configfile path/to/config.yaml
```

### Direct Python Entry

The old Python entrypoint still works:

```bash
python main.py
python main.py --step qc
python main.py --step deseq
python main.py --step gsea
python main.py --step motif
python main.py --step report
```

You can also point it to a YAML config:

```bash
RNASEQ_CONFIG=workflow_config.yaml python main.py
```

## Outputs

Outputs are written under `OUTPUT_DIR`, by default `Final_Analysis_Pipeline/`.

Main directories:

- `00_Validation/`
  contains validated input markers plus sample-group and contrast summaries
- `01_QC/`
  contains original and adjusted PCA outputs, correlation/distance plots, library-size and detected-gene plots, transformed matrices, sample QC metrics, and `QC_Metadata_Associations.tsv`
- `02_DESeq2_Stats/`
  contains one CSV per contrast plus `_contrast_summary.csv`
- `03_Volcano_Plots/`
  contains one labeled volcano plot and one labeled MA plot per contrast
- `04_GSEA/`
  contains one subdirectory per contrast and gene-set library
- `05_Summary/`
  contains `Master_Expression_Table.csv`, `Analysis_Summary.md`, `Report_Index.md`, `Report_Index.html`, annotated heatmaps, DEG count summaries, sample outlier report, QC adjustment comparison, and GSEA summary/clustering outputs

`Report_Index.html` is the recommended delivery entrypoint for browsing the final results.
- `06_Motif/`
  contains optional HOMER results

## Notes On QC: `vst` vs `log1p`

QC now defaults to `vst`.

Why:

- raw counts have strong mean-variance dependence
- simple `log1p` reduces this only partially
- `vst` usually gives more stable PCA and correlation structure for RNA-seq samples

When to switch:

- use `vst` for most datasets
- use `log1p` only if you want a lighter-weight fallback

`rlog` is another variance-stabilizing approach known from DESeq2, but it is slower and not the current default here.

## Notes On LFC Shrinkage

The workflow now supports optional LFC shrinkage.

Why it helps:

- raw LFC estimates are noisier for low-count genes
- shrinkage improves effect-size stability for ranking and visualization
- p-values and adjusted p-values remain based on the original test

Current behavior:

- shrinkage is attempted only when the contrast is directly against that factor's reference level
- non-reference-vs-non-reference contrasts are left unshrunk to avoid misleading coefficient mapping
- when shrinkage is applied, the result table keeps both raw and shrunk LFC columns

## Current Guardrails

The workflow now fails fast on several common problems:

- duplicate sample IDs in metadata
- missing values in the design column
- reference level not present in metadata
- contrasts that reference nonexistent groups
- sample mismatches between metadata and count matrix

## Optional HOMER Setup

If you enable motif analysis, make sure `findMotifs.pl` is available in the environment used by that rule. Depending on your installation, you may also need:

```bash
perl $(which configureHomer.pl) -install mouse
```

or:

```bash
perl $(which configureHomer.pl) -install human
```

## What Was Recently Improved

- workflow execution is now Snakemake-first
- rule environments are managed with `--use-conda`
- formula-based designs and optional contrast tables are supported
- QC uses `vst` by default
- GSEA ranking defaults to `stat`
- LFC shrinkage is supported for reference-based contrasts
- DESeq2 result loading is restricted to configured contrasts
- a per-contrast DEG summary table is generated
- input validation is stricter and fails earlier

## Next Reasonable Improvements

- add optional batch-aware QC summaries and covariate diagnostic plots
- add gene annotation harmonization and ID conversion utilities before GSEA/motif
- add richer final reports and combined QC summary outputs
- split runtime environments further if motif dependencies diverge strongly from the DE/GSEA stack

## What Is Still Missing In A Strong RNA-seq Downstream

This workflow is now solid for standard DE-driven downstream analysis, but a stronger production-grade RNA-seq stack often also includes:

- MA plots, sample-level QC metrics, and covariate association summaries
  these are already included here now
- replicate-aware heatmaps for top variable genes or top DE genes
  still missing
- batch-effect visualizations before and after adjustment
  still missing
- explicit outlier detection reports based on sample distance or Cook's diagnostics
  sample-distance-driven and Cook's-based outlier reporting are now included
- gene annotation harmonization and ID conversion for Ensembl, symbol, Entrez consistency
  partially addressed here by optional gene-version stripping, but not full annotation mapping
- pathway summary tables merged across contrasts
  cross-contrast GSEA summary and term clustering are now included; full enrichment-map style graph summaries are still missing
- publication-grade summary figures such as DEG count bars, upset plots, and contrast-cluster heatmaps
  still missing

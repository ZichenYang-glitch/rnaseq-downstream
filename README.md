# RNA-seq Downstream Analysis Pipeline

A modular, Python-based pipeline for automated RNA-seq downstream analysis. It integrates Quality Control (QC), Differential Expression Analysis (PyDESeq2), Functional Enrichment (GSEA), Transcription Factor Motif Analysis (HOMER), and comprehensive reporting into a single workflow.

## 🚀 Features

*   **Automated Data Loading**: Reads Salmon gene counts and TPMs directly.
*   **Quality Control**: Generates PCA plots and Sample Correlation Heatmaps.
*   **Differential Expression**: Uses `PyDESeq2` to model expression and compute statistics for custom contrasts.
*   **Visualization**: Automatically generates Volcano Plots for every comparison.
*   **Enrichment Analysis**: Runs GSEA Preranked (using `GSEAPy`) for KEGG/GO pathways.
*   **Motif Analysis**: Integrated **HOMER** support to identify enriched transcription factor binding sites in the promoters of DEGs.
*   **Master Reporting**: Aggregates all stats and expression values into a single summary CSV.

## 🛠 Prerequisites

### Installation via Conda (Recommended)

1.  Clone the repository:
    ```bash
    git clone https://github.com/ZichenYang-glitch/rnaseq-downstream.git
    cd rnaseq-downstream
    ```

2.  Create the environment:
    ```bash
    conda env create -f environment.yaml
    ```

3.  Activate the environment:
    ```bash
    conda activate rnaseq-downstream
    ```

### HOMER (Optional)
For Motif analysis, you may need to manually configure HOMER after installation:
```bash
# If you uncommented 'homer' in environment.yaml
perl $(which configureHomer.pl) -install mouse  # or human
```

## 📂 Project Structure

```text
rnaseq-downstream/
├── main.py             # Entry point (Run this!)
├── config.py           # User configuration (Edit this!)
├── README.md           # Documentation
└── modules/            # Analysis Logic
    ├── data.py         # Data loading & cleaning
    ├── deseq.py        # PyDESeq2 wrapper & QC
    ├── enrichment.py   # GSEA analysis
    ├── motif.py        # HOMER Motif analysis wrapper
    └── report.py       # Master table generation
```

## ⚙️ Configuration

Open `config.py` to customize your analysis:

1.  **Input Files**: Set paths to your `counts.tsv`, `tpm.tsv`, and `metadata.txt`.
2.  **Contrasts**: Define which groups to compare.
3.  **Thresholds**: Adjust P-value (`PADJ_THRESH`) and Log2FC (`LOGFC_THRESH`) cutoffs.
4.  **Motif Settings**: Set `HOMER_SPECIES` (e.g., 'mouse' or 'human') and toggle `RUN_MOTIF`.

## ▶️ Usage

### 1. Run Everything (Default)
```bash
python main.py
```

### 2. Run Specific Steps
Use the `--step` argument to run only parts of the pipeline:

*   **QC Only**: `python main.py --step qc`
*   **DESeq2 Only**: `python main.py --step deseq`
*   **GSEA Only**: `python main.py --step gsea`
*   **Motif Only**: `python main.py --step motif`

## 📊 Outputs

Results are saved in `Final_Analysis_Pipeline/`:

*   `01_QC/`: PCA plots and Correlation heatmaps.
*   `02_DESeq2_Stats/`: CSV files for each contrast.
*   `03_Volcano_Plots/`: Volcano plots.
*   `04_GSEA/`: GSEA enrichment plots and tables.
*   `05_Summary/`: `Master_Expression_Table.csv` (Aggregated stats).
*   `06_Motif/`: HOMER output folders (UP/DOWN regulated genes per contrast).

## 📝 Metadata Format

Ensure your `metadata.txt` is tab-separated:

| sample_id | group |
| :--- | :--- |
| SRX...1 | Naive |
| SRX...2 | Th1_24h |
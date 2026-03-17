import csv
import os


def load_contrast_names():
    default_factor = config["DESIGN_FACTOR"]
    path = config.get("CONTRASTS_FILE")
    if path:
        with open(path, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            names = []
            for row in reader:
                factor = row.get("factor") or default_factor
                treatment = row["treatment"]
                control = row["control"]
                name = row.get("name")
                base = f"{treatment}_vs_{control}"
                names.append(name or (base if factor == default_factor else f"{factor}__{base}"))
            return names
    return [f"{treatment}_vs_{control}" for treatment, control in config["CONTRASTS"]]


OUTPUT_DIR = config["OUTPUT_DIR"]
CONTRAST_NAMES = load_contrast_names()
COUNTS_FILE = config["COUNTS_FILE"]
TPM_FILE = config["TPM_FILE"]
METADATA_FILE = config["METADATA_FILE"]
TPM_INPUTS = [TPM_FILE] if os.path.exists(TPM_FILE) else []
CONFIG_INPUTS = [workflow.configfiles[0]]
if config.get("CONTRASTS_FILE"):
    CONFIG_INPUTS.append(config["CONTRASTS_FILE"])
if config.get("UPSTREAM_MANIFEST"):
    CONFIG_INPUTS.append(config["UPSTREAM_MANIFEST"])
VALIDATION_OUTPUTS = [
    f"{OUTPUT_DIR}/00_Validation/validated_inputs.txt",
    f"{OUTPUT_DIR}/00_Validation/group_summary.tsv",
    f"{OUTPUT_DIR}/00_Validation/contrast_summary.tsv",
    f"{OUTPUT_DIR}/00_Validation/upstream_provenance.tsv",
]
QC_OUTPUTS = [
    f"{OUTPUT_DIR}/01_QC/PCA_Plot.png",
    f"{OUTPUT_DIR}/01_QC/Adjusted_PCA_Plot.png",
    f"{OUTPUT_DIR}/01_QC/Sample_Correlation.png",
    f"{OUTPUT_DIR}/01_QC/Sample_Distance.png",
    f"{OUTPUT_DIR}/01_QC/Library_Size.png",
    f"{OUTPUT_DIR}/01_QC/Detected_Genes.png",
    f"{OUTPUT_DIR}/01_QC/QC_Transformed_Counts.tsv",
    f"{OUTPUT_DIR}/01_QC/QC_Adjusted_Transformed_Counts.tsv",
    f"{OUTPUT_DIR}/01_QC/PCA_Coordinates.tsv",
    f"{OUTPUT_DIR}/01_QC/PCA_Explained_Variance.tsv",
    f"{OUTPUT_DIR}/01_QC/Adjusted_PCA_Coordinates.tsv",
    f"{OUTPUT_DIR}/01_QC/Adjusted_PCA_Explained_Variance.tsv",
    f"{OUTPUT_DIR}/01_QC/Sample_Correlation.tsv",
    f"{OUTPUT_DIR}/01_QC/Sample_Distance.tsv",
    f"{OUTPUT_DIR}/01_QC/Sample_QC_Metrics.tsv",
    f"{OUTPUT_DIR}/01_QC/QC_Metadata_Associations.tsv",
]
DESEQ_OUTPUTS = [f"{OUTPUT_DIR}/02_DESeq2_Stats/{name}.csv" for name in CONTRAST_NAMES]
DESEQ_SUMMARY = f"{OUTPUT_DIR}/02_DESeq2_Stats/_contrast_summary.csv"
DESEQ_DIAGNOSTICS = [
    f"{OUTPUT_DIR}/02_DESeq2_Stats/_cooks_matrix.tsv",
    f"{OUTPUT_DIR}/02_DESeq2_Stats/_cooks_gene_report.tsv",
    f"{OUTPUT_DIR}/02_DESeq2_Stats/_cooks_sample_report.tsv",
]
VOLCANO_OUTPUTS = [f"{OUTPUT_DIR}/03_Volcano_Plots/Volcano_{name}.png" for name in CONTRAST_NAMES]
MA_OUTPUTS = [f"{OUTPUT_DIR}/03_Volcano_Plots/MA_{name}.png" for name in CONTRAST_NAMES]
GSEA_OUTPUTS = [
    directory(f"{OUTPUT_DIR}/04_GSEA/{name}/{gene_set}")
    for name in CONTRAST_NAMES
    for gene_set in config["GSEA_GENE_SETS"]
] if config.get("RUN_GSEA", True) else []
GSEA_TARGETS = [
    f"{OUTPUT_DIR}/04_GSEA/{name}/{gene_set}"
    for name in CONTRAST_NAMES
    for gene_set in config["GSEA_GENE_SETS"]
] if config.get("RUN_GSEA", True) else []
REPORT_OUTPUT = f"{OUTPUT_DIR}/05_Summary/Master_Expression_Table.csv"
REPORT_SUMMARY_OUTPUT = f"{OUTPUT_DIR}/05_Summary/Analysis_Summary.md"
REPORT_EXTRA_OUTPUTS = [
    f"{OUTPUT_DIR}/05_Summary/Top_Variable_Genes_Heatmap.png",
    f"{OUTPUT_DIR}/05_Summary/Top_DE_Genes_Heatmap.png",
    f"{OUTPUT_DIR}/05_Summary/Annotated_Top_Variable_Clustermap.png",
    f"{OUTPUT_DIR}/05_Summary/DEG_Counts.tsv",
    f"{OUTPUT_DIR}/05_Summary/DEG_Counts_Barplot.png",
    f"{OUTPUT_DIR}/05_Summary/Sample_Outlier_Report.tsv",
    f"{OUTPUT_DIR}/05_Summary/GSEA_Summary.tsv",
    f"{OUTPUT_DIR}/05_Summary/GSEA_Summary_Dotplot.png",
    f"{OUTPUT_DIR}/05_Summary/GSEA_Term_Contrast_Matrix.tsv",
    f"{OUTPUT_DIR}/05_Summary/GSEA_Term_Clustermap.png",
    f"{OUTPUT_DIR}/05_Summary/QC_Adjustment_Comparison.tsv",
    f"{OUTPUT_DIR}/05_Summary/Report_Index.md",
    f"{OUTPUT_DIR}/05_Summary/Report_Index.html",
]
REPORT_CONTRAST_PAGES = [f"{OUTPUT_DIR}/05_Summary/contrast_{name}.html" for name in CONTRAST_NAMES]
MOTIF_OUTPUTS = [
    directory(f"{OUTPUT_DIR}/06_Motif/{name}")
    for name in CONTRAST_NAMES
] if config.get("RUN_MOTIF", False) else []
MOTIF_TARGETS = [
    f"{OUTPUT_DIR}/06_Motif/{name}"
    for name in CONTRAST_NAMES
] if config.get("RUN_MOTIF", False) else []


rule all:
    input:
        VALIDATION_OUTPUTS,
        QC_OUTPUTS,
        DESEQ_OUTPUTS,
        DESEQ_SUMMARY,
        DESEQ_DIAGNOSTICS,
        VOLCANO_OUTPUTS,
        MA_OUTPUTS,
        REPORT_OUTPUT,
        REPORT_SUMMARY_OUTPUT,
        REPORT_EXTRA_OUTPUTS,
        REPORT_CONTRAST_PAGES,
        GSEA_TARGETS,
        MOTIF_TARGETS,


rule validate:
    input:
        COUNTS_FILE,
        METADATA_FILE,
        CONFIG_INPUTS
    output:
        VALIDATION_OUTPUTS
    conda:
        "../../environment.yaml"
    shell:
        "MPLCONFIGDIR=/tmp/mplconfig RNASEQ_CONFIG={workflow.configfiles[0]} python scripts/run_validate.py"


rule qc:
    input:
        COUNTS_FILE,
        METADATA_FILE,
        CONFIG_INPUTS,
        VALIDATION_OUTPUTS
    output:
        QC_OUTPUTS
    conda:
        "../../environment.yaml"
    shell:
        "MPLCONFIGDIR=/tmp/mplconfig RNASEQ_CONFIG={workflow.configfiles[0]} python scripts/run_qc.py"


rule deseq:
    input:
        COUNTS_FILE,
        METADATA_FILE,
        CONFIG_INPUTS,
        VALIDATION_OUTPUTS,
        QC_OUTPUTS
    output:
        DESEQ_OUTPUTS,
        DESEQ_SUMMARY,
        DESEQ_DIAGNOSTICS,
        VOLCANO_OUTPUTS,
        MA_OUTPUTS
    conda:
        "../../environment.yaml"
    shell:
        "MPLCONFIGDIR=/tmp/mplconfig RNASEQ_CONFIG={workflow.configfiles[0]} python scripts/run_deseq.py"


rule gsea:
    input:
        DESEQ_OUTPUTS,
        CONFIG_INPUTS
    output:
        GSEA_OUTPUTS
    conda:
        "../../environment.yaml"
    shell:
        "MPLCONFIGDIR=/tmp/mplconfig RNASEQ_CONFIG={workflow.configfiles[0]} python scripts/run_gsea.py"


rule report:
    input:
        COUNTS_FILE,
        METADATA_FILE,
        TPM_INPUTS,
        QC_OUTPUTS,
        DESEQ_OUTPUTS,
        GSEA_TARGETS,
        VALIDATION_OUTPUTS,
        CONFIG_INPUTS
    output:
        REPORT_OUTPUT,
        REPORT_SUMMARY_OUTPUT,
        REPORT_EXTRA_OUTPUTS,
        REPORT_CONTRAST_PAGES
    conda:
        "../../environment.yaml"
    shell:
        "MPLCONFIGDIR=/tmp/mplconfig RNASEQ_CONFIG={workflow.configfiles[0]} python scripts/run_report.py"


rule motif:
    input:
        DESEQ_OUTPUTS,
        CONFIG_INPUTS
    output:
        MOTIF_OUTPUTS
    conda:
        "../../environment.yaml"
    shell:
        "MPLCONFIGDIR=/tmp/mplconfig RNASEQ_CONFIG={workflow.configfiles[0]} python scripts/run_motif.py"

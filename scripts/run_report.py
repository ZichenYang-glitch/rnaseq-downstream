#!/usr/bin/env python3
import os

import pandas as pd

import config as cfg
from modules import data, report


def main():
    meta = data.load_metadata(cfg.METADATA_FILE, cfg.DESIGN_FACTOR)
    contrasts = data.load_contrasts(cfg.DESIGN_FACTOR, cfg.CONTRASTS, cfg.CONTRASTS_FILE)
    data.validate_analysis_inputs(meta, cfg.DESIGN_FACTOR, cfg.REFERENCE_LEVELS, contrasts)
    meta = data.prepare_metadata(meta, cfg.REFERENCE_LEVELS, cfg.CONTINUOUS_FACTORS)
    upstream_manifest = data.load_upstream_manifest(cfg.UPSTREAM_MANIFEST)
    upstream_provenance = data.build_upstream_provenance(vars(cfg), upstream_manifest)
    annotation_df = data.load_annotation_table(
        cfg.ANNOTATION_FILE,
        gene_id_col=cfg.ANNOTATION_GENE_ID_COL,
        gene_name_col=cfg.ANNOTATION_GENE_NAME_COL,
        strip_gene_version=cfg.STRIP_GENE_VERSION,
    )
    counts_T = data.load_counts(
        cfg.COUNTS_FILE,
        meta.index,
        cfg.MIN_COUNTS,
        strip_gene_version=cfg.STRIP_GENE_VERSION,
    )

    results = {}
    deseq_dir = os.path.join(cfg.OUTPUT_DIR, '02_DESeq2_Stats')
    for contrast in contrasts:
        name = data.contrast_name(contrast, cfg.DESIGN_FACTOR)
        path = os.path.join(deseq_dir, f"{name}.csv")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing DESeq2 result: {path}")
        results[name] = pd.read_csv(path, index_col=0)

    report.create_master_table(
        meta,
        cfg.DESIGN_FACTOR,
        results,
        cfg.TPM_FILE,
        counts_T,
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
        strip_gene_version=cfg.STRIP_GENE_VERSION,
        annotation_df=annotation_df,
    )
    report.create_analysis_summary(
        meta,
        cfg.DESIGN_FACTOR,
        results,
        os.path.join(cfg.OUTPUT_DIR, '01_QC'),
        os.path.join(cfg.OUTPUT_DIR, '00_Validation'),
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
        upstream_provenance=upstream_provenance,
    )
    report.create_heatmap_summaries(
        meta,
        cfg.DESIGN_FACTOR,
        results,
        os.path.join(cfg.OUTPUT_DIR, '01_QC'),
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
        top_variable_genes=cfg.TOP_VARIABLE_GENES,
        top_de_genes_per_contrast=cfg.TOP_DE_GENES_PER_CONTRAST,
        annotation_factors=cfg.QC_ANNOTATION_FACTORS,
    )
    report.create_deg_summary_plot(
        results,
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
        padj_thresh=cfg.PADJ_THRESH,
        logfc_thresh=cfg.LOGFC_THRESH,
    )
    report.create_sample_outlier_report(
        os.path.join(cfg.OUTPUT_DIR, '01_QC'),
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
    )
    report.create_gsea_summary(
        os.path.join(cfg.OUTPUT_DIR, '04_GSEA'),
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
        fdr_thresh=cfg.GSEA_SUMMARY_FDR,
        top_terms=cfg.GSEA_SUMMARY_TOP_TERMS,
    )
    report.create_qc_adjustment_comparison(
        os.path.join(cfg.OUTPUT_DIR, '01_QC'),
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
    )
    report.create_report_index(
        cfg.OUTPUT_DIR,
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
        upstream_provenance=upstream_provenance,
    )
    report.create_html_report_index(
        cfg.OUTPUT_DIR,
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
        upstream_provenance=upstream_provenance,
    )
    report.create_contrast_report_pages(
        cfg.OUTPUT_DIR,
        os.path.join(cfg.OUTPUT_DIR, '05_Summary'),
        results,
    )


if __name__ == "__main__":
    main()

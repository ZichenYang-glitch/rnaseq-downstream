#!/usr/bin/env python3
import argparse
import sys
import os
import matplotlib
matplotlib.use('Agg') # Set non-interactive backend
import config as cfg

# Import modules
from modules import data, deseq, enrichment, report, motif

def main():
    parser = argparse.ArgumentParser(description="RNA-seq Downstream Analysis Pipeline")
    parser.add_argument('--step', type=str, choices=['all', 'qc', 'deseq', 'gsea', 'motif', 'report'], 
                        default='all', help="Analysis step to run.")
    args = parser.parse_args()
    
    print("=== RNA-seq Pipeline Started ===")
    
    # 1. Load Data (Always required)
    print("[1] Loading Data...")
    meta = data.load_metadata(cfg.METADATA_FILE, cfg.DESIGN_FACTOR)
    contrasts = data.load_contrasts(cfg.DESIGN_FACTOR, cfg.CONTRASTS, cfg.CONTRASTS_FILE)
    data.validate_analysis_inputs(meta, cfg.DESIGN_FACTOR, cfg.REFERENCE_LEVELS, contrasts)
    meta = data.prepare_metadata(meta, cfg.REFERENCE_LEVELS, cfg.CONTINUOUS_FACTORS)
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
    print(f"    Loaded {counts_T.shape[0]} samples, {counts_T.shape[1]} genes.")
    
    # Paths
    dirs = {
        'qc': os.path.join(cfg.OUTPUT_DIR, '01_QC'),
        'deseq': os.path.join(cfg.OUTPUT_DIR, '02_DESeq2_Stats'),
        'plots': os.path.join(cfg.OUTPUT_DIR, '03_Volcano_Plots'),
        'gsea': os.path.join(cfg.OUTPUT_DIR, '04_GSEA'),
        'motif': os.path.join(cfg.OUTPUT_DIR, '06_Motif'),
        'report': os.path.join(cfg.OUTPUT_DIR, '05_Summary')
    }
    
    # 2. QC
    if args.step in ['all', 'qc']:
        print("[2] Running QC...")
        deseq.run_qc(
            counts_T,
            meta,
            cfg.DESIGN_FACTOR,
            dirs['qc'],
            design=cfg.DESIGN,
            continuous_factors=cfg.CONTINUOUS_FACTORS,
            transform=cfg.QC_TRANSFORM,
            adjust_factors=cfg.QC_ADJUST_FACTORS,
            use_design=cfg.VST_USE_DESIGN,
            label_samples=cfg.QC_LABEL_SAMPLES,
            n_cpus=cfg.N_CPUS,
        )
        if args.step == 'qc': sys.exit(0)
        
    # 3. DESeq2
    results = {}
    if args.step in ['all', 'deseq', 'gsea', 'motif', 'report']:
        # Check if we can load existing results or need to run
        if args.step == 'deseq' or args.step == 'all' or not os.path.exists(dirs['deseq']):
            print("[3] Running DESeq2...")
            # Ensure metadata matches counts samples exactly and in the same order
            meta_subset = meta.loc[counts_T.index]
            dds = deseq.fit_model(
                counts_T,
                meta_subset,
                cfg.DESIGN,
                cfg.CONTINUOUS_FACTORS,
                n_cpus=cfg.N_CPUS,
            )
            results = deseq.run_contrasts(
                dds,
                contrasts,
                cfg.DESIGN_FACTOR,
                dirs['deseq'],
                n_cpus=cfg.N_CPUS,
                reference_levels=cfg.REFERENCE_LEVELS,
                shrink_lfc=cfg.SHRINK_LFC,
                shrink_adapt=cfg.SHRINK_LFC_ADAPT,
            )
            deseq.write_cooks_diagnostics(dds, dirs['deseq'])
            deseq.write_contrast_summary(results, dirs['deseq'], cfg.PADJ_THRESH, cfg.LOGFC_THRESH)
            deseq.plot_volcano(results, dirs['plots'], cfg.PADJ_THRESH, cfg.LOGFC_THRESH, cfg.TOP_LABEL_GENES)
            deseq.plot_ma(results, dirs['plots'], cfg.PADJ_THRESH, cfg.LOGFC_THRESH, cfg.TOP_LABEL_GENES)
        else:
            print("[3] Loading existing DESeq2 results...")
            import pandas as pd
            for contrast in contrasts:
                name = data.contrast_name(contrast, cfg.DESIGN_FACTOR)
                path = os.path.join(dirs['deseq'], f"{name}.csv")
                if not os.path.exists(path):
                    raise FileNotFoundError(
                        f"Expected DESeq2 result not found for contrast '{name}': {path}"
                    )
                results[name] = pd.read_csv(path, index_col=0)
                
        if args.step == 'deseq': sys.exit(0)
        
    # 4. GSEA
    if args.step in ['all', 'gsea'] and cfg.RUN_GSEA:
        print("[4] Running GSEA...")
        enrichment.run_gsea(
            results,
            cfg.GSEA_GENE_SETS,
            dirs['gsea'],
            cfg.GSEA_PERMUTATIONS,
            cfg.GSEA_RANK_METRIC,
        )
        
    # 5. Motif Analysis
    if args.step in ['all', 'motif'] and cfg.RUN_MOTIF:
        print("[5] Running Motif Analysis (HOMER)...")
        motif.run_motif_analysis(results, dirs['motif'], cfg.HOMER_SPECIES, cfg.PADJ_THRESH, cfg.LOGFC_THRESH)

    # 6. Report
    if args.step in ['all', 'report']:
        print("[6] Generating Report...")
        report.create_master_table(
            meta,
            cfg.DESIGN_FACTOR,
            results,
            cfg.TPM_FILE,
            counts_T,
            dirs['report'],
            strip_gene_version=cfg.STRIP_GENE_VERSION,
            annotation_df=annotation_df,
        )
        report.create_analysis_summary(
            meta,
            cfg.DESIGN_FACTOR,
            results,
            dirs['qc'],
            os.path.join(cfg.OUTPUT_DIR, '00_Validation'),
            dirs['report'],
        )
        report.create_heatmap_summaries(
            meta,
            cfg.DESIGN_FACTOR,
            results,
            dirs['qc'],
            dirs['report'],
            top_variable_genes=cfg.TOP_VARIABLE_GENES,
            top_de_genes_per_contrast=cfg.TOP_DE_GENES_PER_CONTRAST,
            annotation_factors=cfg.QC_ANNOTATION_FACTORS,
        )
        report.create_deg_summary_plot(
            results,
            dirs['report'],
            padj_thresh=cfg.PADJ_THRESH,
            logfc_thresh=cfg.LOGFC_THRESH,
        )
        report.create_sample_outlier_report(dirs['qc'], dirs['report'])
        report.create_gsea_summary(
            dirs['gsea'],
            dirs['report'],
            fdr_thresh=cfg.GSEA_SUMMARY_FDR,
            top_terms=cfg.GSEA_SUMMARY_TOP_TERMS,
        )
        report.create_qc_adjustment_comparison(dirs['qc'], dirs['report'])
        report.create_report_index(cfg.OUTPUT_DIR, dirs['report'])
        report.create_html_report_index(cfg.OUTPUT_DIR, dirs['report'])
        report.create_contrast_report_pages(cfg.OUTPUT_DIR, dirs['report'], results)
        
    print("\n=== Pipeline Complete ===")
    print(f"Results located in: {os.path.abspath(cfg.OUTPUT_DIR)}")

if __name__ == "__main__":
    main()

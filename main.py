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
    counts_T = data.load_counts(cfg.COUNTS_FILE, meta.index, cfg.MIN_COUNTS)
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
        deseq.run_qc(counts_T, meta, cfg.DESIGN_FACTOR, dirs['qc'])
        if args.step == 'qc': sys.exit(0)
        
    # 3. DESeq2
    results = {}
    if args.step in ['all', 'deseq', 'gsea', 'report']:
        # Check if we can load existing results or need to run
        ran_deseq = False
        if args.step == 'deseq' or args.step == 'all' or not os.path.exists(dirs['deseq']):
            print("[3] Running DESeq2...")
            # Ensure metadata matches counts samples exactly and in the same order
            meta_subset = meta.loc[counts_T.index]
            dds = deseq.fit_model(counts_T, meta_subset, cfg.DESIGN_FACTOR, cfg.REFERENCE_LEVEL)
            results = deseq.run_contrasts(dds, cfg.CONTRASTS, cfg.DESIGN_FACTOR, dirs['deseq'])
            deseq.plot_volcano(results, dirs['plots'], cfg.PADJ_THRESH, cfg.LOGFC_THRESH)
            ran_deseq = True
        else:
            print("[3] Loading existing DESeq2 results...")
            # Load from CSVs
            import pandas as pd
            import glob
            files = glob.glob(os.path.join(dirs['deseq'], "*.csv"))
            for f in files:
                name = os.path.basename(f).replace('.csv', '')
                results[name] = pd.read_csv(f, index_col=0)
                
        if args.step == 'deseq': sys.exit(0)
        
    # 4. GSEA
    if args.step in ['all', 'gsea'] and cfg.RUN_GSEA:
        print("[4] Running GSEA...")
        enrichment.run_gsea(results, cfg.GSEA_GENE_SETS, dirs['gsea'])
        
    # 5. Motif Analysis
    if args.step in ['all', 'motif'] and cfg.RUN_MOTIF:
        print("[5] Running Motif Analysis (HOMER)...")
        motif.run_motif_analysis(results, dirs['motif'], cfg.HOMER_SPECIES, cfg.PADJ_THRESH, cfg.LOGFC_THRESH)

    # 6. Report
    if args.step in ['all', 'report']:
        print("[6] Generating Report...")
        report.create_master_table(meta, cfg.DESIGN_FACTOR, results, cfg.TPM_FILE, counts_T, dirs['report'])
        
    print("\n=== Pipeline Complete ===")
    print(f"Results located in: {os.path.abspath(cfg.OUTPUT_DIR)}")

if __name__ == "__main__":
    main()

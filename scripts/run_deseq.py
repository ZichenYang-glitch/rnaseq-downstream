#!/usr/bin/env python3
import os

import config as cfg
from modules import data, deseq


def main():
    meta = data.load_metadata(cfg.METADATA_FILE, cfg.DESIGN_FACTOR)
    contrasts = data.load_contrasts(cfg.DESIGN_FACTOR, cfg.CONTRASTS, cfg.CONTRASTS_FILE)
    data.validate_analysis_inputs(meta, cfg.DESIGN_FACTOR, cfg.REFERENCE_LEVELS, contrasts)
    meta = data.prepare_metadata(meta, cfg.REFERENCE_LEVELS, cfg.CONTINUOUS_FACTORS)
    counts_T = data.load_counts(
        cfg.COUNTS_FILE,
        meta.index,
        cfg.MIN_COUNTS,
        strip_gene_version=cfg.STRIP_GENE_VERSION,
    )
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
        os.path.join(cfg.OUTPUT_DIR, '02_DESeq2_Stats'),
        n_cpus=cfg.N_CPUS,
        reference_levels=cfg.REFERENCE_LEVELS,
        shrink_lfc=cfg.SHRINK_LFC,
        shrink_adapt=cfg.SHRINK_LFC_ADAPT,
    )
    deseq.write_cooks_diagnostics(dds, os.path.join(cfg.OUTPUT_DIR, '02_DESeq2_Stats'))
    deseq.write_contrast_summary(
        results,
        os.path.join(cfg.OUTPUT_DIR, '02_DESeq2_Stats'),
        cfg.PADJ_THRESH,
        cfg.LOGFC_THRESH,
    )
    deseq.plot_volcano(
        results,
        os.path.join(cfg.OUTPUT_DIR, '03_Volcano_Plots'),
        cfg.PADJ_THRESH,
        cfg.LOGFC_THRESH,
        cfg.TOP_LABEL_GENES,
    )
    deseq.plot_ma(
        results,
        os.path.join(cfg.OUTPUT_DIR, '03_Volcano_Plots'),
        cfg.PADJ_THRESH,
        cfg.LOGFC_THRESH,
        cfg.TOP_LABEL_GENES,
    )


if __name__ == "__main__":
    main()

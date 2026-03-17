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
    out_dir = os.path.join(cfg.OUTPUT_DIR, '01_QC')
    deseq.run_qc(
        counts_T,
        meta,
        cfg.DESIGN_FACTOR,
        out_dir,
        design=cfg.DESIGN,
        continuous_factors=cfg.CONTINUOUS_FACTORS,
        transform=cfg.QC_TRANSFORM,
        adjust_factors=cfg.QC_ADJUST_FACTORS,
        use_design=cfg.VST_USE_DESIGN,
        label_samples=cfg.QC_LABEL_SAMPLES,
        n_cpus=cfg.N_CPUS,
    )


if __name__ == "__main__":
    main()

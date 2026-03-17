#!/usr/bin/env python3
import os

import pandas as pd

import config as cfg
from modules import data, motif


def main():
    contrasts = data.load_contrasts(cfg.DESIGN_FACTOR, cfg.CONTRASTS, cfg.CONTRASTS_FILE)
    results = {}
    deseq_dir = os.path.join(cfg.OUTPUT_DIR, '02_DESeq2_Stats')
    for contrast in contrasts:
        name = data.contrast_name(contrast, cfg.DESIGN_FACTOR)
        path = os.path.join(deseq_dir, f"{name}.csv")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing DESeq2 result: {path}")
        results[name] = pd.read_csv(path, index_col=0)

    motif.run_motif_analysis(
        results,
        os.path.join(cfg.OUTPUT_DIR, '06_Motif'),
        cfg.HOMER_SPECIES,
        cfg.PADJ_THRESH,
        cfg.LOGFC_THRESH,
    )


if __name__ == "__main__":
    main()

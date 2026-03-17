#!/usr/bin/env python3
import os

import pandas as pd

import config as cfg
from modules import data


def main():
    meta = data.load_metadata(cfg.METADATA_FILE, cfg.DESIGN_FACTOR)
    contrasts = data.load_contrasts(cfg.DESIGN_FACTOR, cfg.CONTRASTS, cfg.CONTRASTS_FILE)
    data.validate_analysis_inputs(meta, cfg.DESIGN_FACTOR, cfg.REFERENCE_LEVELS, contrasts)
    meta = data.prepare_metadata(meta, cfg.REFERENCE_LEVELS, cfg.CONTINUOUS_FACTORS)
    upstream_manifest = data.load_upstream_manifest(cfg.UPSTREAM_MANIFEST)
    upstream_provenance = data.build_upstream_provenance(vars(cfg), upstream_manifest)
    counts_T = data.load_counts(
        cfg.COUNTS_FILE,
        meta.index,
        cfg.MIN_COUNTS,
        strip_gene_version=cfg.STRIP_GENE_VERSION,
    )

    out_dir = os.path.join(cfg.OUTPUT_DIR, '00_Validation')
    os.makedirs(out_dir, exist_ok=True)

    group_summary = (
        meta[cfg.DESIGN_FACTOR]
        .astype(str)
        .value_counts()
        .rename_axis(cfg.DESIGN_FACTOR)
        .reset_index(name='n_samples')
        .sort_values(cfg.DESIGN_FACTOR)
    )
    group_summary.to_csv(os.path.join(out_dir, 'group_summary.tsv'), sep='\t', index=False)

    contrast_summary = pd.DataFrame([
        {
            'name': data.contrast_name(contrast, cfg.DESIGN_FACTOR),
            'factor': contrast['factor'],
            'treatment': contrast['treatment'],
            'control': contrast['control'],
        }
        for contrast in contrasts
    ])
    contrast_summary.to_csv(os.path.join(out_dir, 'contrast_summary.tsv'), sep='\t', index=False)

    upstream_rows = [{'field': key, 'value': value} for key, value in upstream_provenance.items()]
    if not upstream_rows:
        upstream_rows = [{'field': 'upstream_pipeline_name', 'value': 'not_provided'}]
    pd.DataFrame(upstream_rows).to_csv(
        os.path.join(out_dir, 'upstream_provenance.tsv'),
        sep='\t',
        index=False,
    )

    with open(os.path.join(out_dir, 'validated_inputs.txt'), 'w', encoding='utf-8') as handle:
        handle.write(f"samples\t{counts_T.shape[0]}\n")
        handle.write(f"genes\t{counts_T.shape[1]}\n")
        handle.write(f"design\t{cfg.DESIGN}\n")
        handle.write(f"design_factor\t{cfg.DESIGN_FACTOR}\n")
        handle.write(f"n_contrasts\t{len(contrasts)}\n")


if __name__ == "__main__":
    main()

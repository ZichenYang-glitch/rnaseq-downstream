import os

import yaml

# ================= CONFIGURATION =================

DEFAULTS = {
    'COUNTS_FILE': 'salmon.merged.gene_counts.tsv',
    'TPM_FILE': 'salmon.merged.gene_tpm.tsv',
    'METADATA_FILE': 'metadata.txt',
    'OUTPUT_DIR': 'Final_Analysis_Pipeline',
    'UPSTREAM_MANIFEST': None,
    'UPSTREAM_PIPELINE_NAME': 'Nextflow RNA-seq pipeline',
    'UPSTREAM_PIPELINE_VERSION': None,
    'UPSTREAM_PIPELINE_URL': None,
    'REFERENCE_GENOME': None,
    'ANNOTATION_RELEASE': None,
    'QUANTIFICATION_METHOD': None,
    'COUNT_MATRIX_TYPE': 'gene_counts',
    'TPM_MATRIX_TYPE': 'gene_tpm',
    'STRIP_GENE_VERSION': True,
    'ANNOTATION_FILE': None,
    'ANNOTATION_GENE_ID_COL': 'gene_id',
    'ANNOTATION_GENE_NAME_COL': 'gene_name',
    'DESIGN': None,
    'DESIGN_FACTOR': 'group',
    'REFERENCE_LEVEL': 'M0',
    'REFERENCE_LEVELS': {},
    'CONTINUOUS_FACTORS': [],
    'CONTRASTS_FILE': None,
    'CONTRASTS': [
        ('M1', 'M0'),
        ('M2', 'M0'),
        ('M1', 'M2'),
    ],
    'MIN_COUNTS': 10,
    'PADJ_THRESH': 0.05,
    'LOGFC_THRESH': 1.0,
    'TOP_VARIABLE_GENES': 50,
    'TOP_DE_GENES_PER_CONTRAST': 30,
    'TOP_LABEL_GENES': 10,
    'GSEA_SUMMARY_FDR': 0.25,
    'GSEA_SUMMARY_TOP_TERMS': 20,
    'QC_TRANSFORM': 'vst',
    'QC_ADJUST_FACTORS': [],
    'QC_ANNOTATION_FACTORS': [],
    'VST_USE_DESIGN': False,
    'QC_LABEL_SAMPLES': False,
    'RUN_GSEA': True,
    'GSEA_GENE_SETS': ['KEGG_2019_Mouse', 'GO_Biological_Process_2021'],
    'GSEA_PERMUTATIONS': 1000,
    'GSEA_RANK_METRIC': 'stat',
    'SHRINK_LFC': True,
    'SHRINK_LFC_ADAPT': True,
    'RUN_MOTIF': False,
    'HOMER_SPECIES': 'mouse',
    'N_CPUS': 4,
}


def _load_overrides():
    cfg_path = os.environ.get('RNASEQ_CONFIG')
    if not cfg_path:
        return {}

    with open(cfg_path, 'r', encoding='utf-8') as handle:
        loaded = yaml.safe_load(handle) or {}

    if not isinstance(loaded, dict):
        raise ValueError(f"Config file must define a mapping: {cfg_path}")

    overrides = dict(loaded)
    if 'CONTRASTS' in overrides:
        overrides['CONTRASTS'] = [tuple(item) for item in overrides['CONTRASTS']]
    return overrides


_CONFIG = DEFAULTS | _load_overrides()
if not _CONFIG.get('DESIGN'):
    _CONFIG['DESIGN'] = f"~ {_CONFIG['DESIGN_FACTOR']}"

_reference_levels = dict(_CONFIG.get('REFERENCE_LEVELS', {}))
if _CONFIG.get('DESIGN_FACTOR') and _CONFIG.get('REFERENCE_LEVEL') is not None:
    _reference_levels.setdefault(_CONFIG['DESIGN_FACTOR'], _CONFIG['REFERENCE_LEVEL'])
_CONFIG['REFERENCE_LEVELS'] = _reference_levels

globals().update(_CONFIG)

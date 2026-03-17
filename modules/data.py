import pandas as pd
import sys
import os
import json

import yaml


def normalize_gene_ids(index, strip_gene_version=True):
    values = pd.Index(index).astype(str)
    if strip_gene_version:
        return values.str.replace(r'\.\d+$', '', regex=True)
    return values


def load_annotation_table(
    path,
    gene_id_col='gene_id',
    gene_name_col='gene_name',
    strip_gene_version=True,
):
    if not path:
        return None

    df = pd.read_csv(path, sep=None, engine='python')
    if gene_id_col not in df.columns:
        raise ValueError(f"Annotation file missing gene ID column '{gene_id_col}'.")
    if gene_name_col not in df.columns:
        raise ValueError(f"Annotation file missing gene name column '{gene_name_col}'.")

    ann = df[[gene_id_col, gene_name_col]].copy()
    ann[gene_id_col] = normalize_gene_ids(ann[gene_id_col], strip_gene_version=strip_gene_version)
    ann = ann.dropna().drop_duplicates(subset=[gene_id_col])
    ann.columns = ['gene_id', 'gene_name']
    return ann.set_index('gene_id')


def load_upstream_manifest(path):
    if not path:
        return {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"Upstream manifest not found: {path}")

    suffix = os.path.splitext(path)[1].lower()
    if suffix in {'.yaml', '.yml'}:
        with open(path, 'r', encoding='utf-8') as handle:
            data = yaml.safe_load(handle) or {}
    elif suffix == '.json':
        with open(path, 'r', encoding='utf-8') as handle:
            data = json.load(handle) or {}
    else:
        df = pd.read_csv(path, sep=None, engine='python')
        if df.shape[1] < 2:
            raise ValueError(
                "Tabular upstream manifest must contain at least two columns: key and value."
            )
        data = dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(str)))

    if not isinstance(data, dict):
        raise ValueError(f"Upstream manifest must define a mapping: {path}")
    return data


def build_upstream_provenance(config_values, manifest=None):
    manifest = manifest or {}
    provenance = {
        'upstream_manifest': config_values.get('UPSTREAM_MANIFEST'),
        'upstream_pipeline_name': config_values.get('UPSTREAM_PIPELINE_NAME'),
        'upstream_pipeline_version': config_values.get('UPSTREAM_PIPELINE_VERSION'),
        'upstream_pipeline_url': config_values.get('UPSTREAM_PIPELINE_URL'),
        'reference_genome': config_values.get('REFERENCE_GENOME'),
        'annotation_release': config_values.get('ANNOTATION_RELEASE'),
        'quantification_method': config_values.get('QUANTIFICATION_METHOD'),
        'count_matrix_type': config_values.get('COUNT_MATRIX_TYPE'),
        'tpm_matrix_type': config_values.get('TPM_MATRIX_TYPE'),
    }

    for key, value in manifest.items():
        provenance.setdefault(str(key), value)
        if provenance.get(str(key)) in [None, '']:
            provenance[str(key)] = value

    return {key: value for key, value in provenance.items() if value not in [None, '']}


def load_expression_matrix(path, metadata_samples=None, strip_gene_version=True, aggregate='sum'):
    df = pd.read_csv(path, sep='\t')
    idx = 'gene_name' if 'gene_name' in df.columns else 'gene_id'
    if idx not in df.columns:
        idx = df.columns[0]

    df[idx] = normalize_gene_ids(df[idx], strip_gene_version=strip_gene_version)
    grouped = getattr(df.groupby(idx), aggregate)(numeric_only=True)

    if metadata_samples is None:
        return grouped

    common = [s for s in metadata_samples if s in grouped.columns]
    if not common:
        raise ValueError("No matching samples between metadata and expression matrix!")

    missing = [s for s in metadata_samples if s not in grouped.columns]
    extra = [s for s in grouped.columns if s not in metadata_samples]
    if missing:
        print(f"[Warning] Samples in metadata but not matrix: {', '.join(missing)}")
    if extra:
        print(f"[Warning] Samples in matrix but not metadata: {', '.join(extra)}")

    return grouped[common]


def contrast_name(contrast, default_factor):
    factor = contrast['factor']
    base = f"{contrast['treatment']}_vs_{contrast['control']}"
    if factor == default_factor:
        return contrast.get('name') or base
    return contrast.get('name') or f"{factor}__{base}"


def load_contrasts(default_factor, inline_contrasts=None, contrasts_file=None):
    if contrasts_file:
        df = pd.read_csv(contrasts_file, sep=None, engine='python')
        required = {'treatment', 'control'}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"Contrasts file must contain columns: {', '.join(sorted(required))}. "
                f"Missing: {', '.join(sorted(missing))}"
            )
        records = []
        for row in df.to_dict(orient='records'):
            factor = row.get('factor') or default_factor
            records.append({
                'factor': str(factor),
                'treatment': str(row['treatment']),
                'control': str(row['control']),
                'name': str(row['name']) if row.get('name') else None,
            })
        return records

    records = []
    for item in inline_contrasts or []:
        if isinstance(item, dict):
            records.append({
                'factor': str(item.get('factor', default_factor)),
                'treatment': str(item['treatment']),
                'control': str(item['control']),
                'name': str(item['name']) if item.get('name') else None,
            })
        else:
            treatment, control = item
            records.append({
                'factor': str(default_factor),
                'treatment': str(treatment),
                'control': str(control),
                'name': None,
            })
    return records


def prepare_metadata(metadata, reference_levels, continuous_factors=None):
    df = metadata.copy()
    for factor, ref in (reference_levels or {}).items():
        if factor not in df.columns:
            raise ValueError(f"Reference level provided for unknown metadata column '{factor}'.")
        levels = [str(value) for value in df[factor].dropna().astype(str).unique()]
        if str(ref) not in levels:
            raise ValueError(
                f"Reference level '{ref}' not found in metadata column '{factor}'. "
                f"Available groups: {', '.join(sorted(levels))}"
            )
        ordered_levels = [str(ref)] + [level for level in sorted(levels) if level != str(ref)]
        df[factor] = pd.Categorical(df[factor].astype(str), categories=ordered_levels, ordered=True)

    for factor in continuous_factors or []:
        if factor not in df.columns:
            raise ValueError(f"Continuous factor '{factor}' not found in metadata.")
        df[factor] = pd.to_numeric(df[factor], errors='raise')

    return df


def load_metadata(path, design_col):
    """Loads and validates metadata."""
    try:
        df = pd.read_csv(path, sep=r'\s+', engine='python')
        # Handle index
        if 'sample_id' in df.columns:
            df = df.set_index('sample_id')
        else:
            df = df.set_index(df.columns[0])
            
        if design_col not in df.columns:
            raise ValueError(f"Design column '{design_col}' not found in metadata.")
        if df.index.has_duplicates:
            dup_ids = df.index[df.index.duplicated()].unique().tolist()
            raise ValueError(f"Duplicate sample IDs found in metadata: {', '.join(map(str, dup_ids))}")
        if df[design_col].isna().any():
            raise ValueError(f"Missing values found in design column '{design_col}'.")
            
        return df
    except Exception as e:
        print(f"[Error] Loading metadata: {e}")
        sys.exit(1)

def load_counts(path, metadata_samples, min_counts=10, strip_gene_version=True):
    """Loads counts, aligns with metadata, and filters."""
    try:
        df = load_expression_matrix(
            path,
            metadata_samples,
            strip_gene_version=strip_gene_version,
            aggregate='sum',
        )
        df = df.fillna(0).round().astype(int)
        
        # Filter
        df = df[df.sum(axis=1) >= min_counts]
        
        return df.T  # Return samples x genes
    except Exception as e:
        print(f"[Error] Loading counts: {e}")
        sys.exit(1)


def validate_analysis_inputs(metadata, design_col, reference_levels, contrasts):
    """Fails fast on invalid design settings and contrasts."""
    if design_col not in metadata.columns:
        raise ValueError(f"Design factor '{design_col}' not found in metadata.")

    for factor, ref in (reference_levels or {}).items():
        groups = set(metadata[factor].astype(str).unique())
        if str(ref) not in groups:
            raise ValueError(
                f"Reference level '{ref}' not found in metadata column '{factor}'. "
                f"Available groups: {', '.join(sorted(groups))}"
            )

    invalid = []
    for contrast in contrasts:
        factor = contrast['factor']
        if factor not in metadata.columns:
            invalid.append(contrast_name(contrast, design_col))
            continue
        groups = set(metadata[factor].astype(str).unique())
        if str(contrast['treatment']) not in groups or str(contrast['control']) not in groups:
            invalid.append(contrast_name(contrast, design_col))

    if invalid:
        raise ValueError(
            f"Contrasts contain factors or groups absent from metadata: {', '.join(invalid)}"
        )

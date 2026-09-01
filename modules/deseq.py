import hashlib
import json
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.metrics import pairwise_distances, silhouette_score

from modules import data
from rnaseq_downstream.errors import QCValidationError
from rnaseq_downstream.qc_math import (
    PCA_METHOD_ID,
    FactorSpec,
    build_joint_design,
    centered_unscaled_pca,
    metadata_reorder_indices,
    remove_nuisance_effects,
    resolve_adjustment_biology_factors,
    select_top_variable_features,
    validate_factor_names,
)


_PCA_METHOD_FILE = 'QC_PCA_Method.json'
_PCA_SELECTED_GENES_FILE = 'QC_PCA_Selected_Genes.tsv'
_CORRECTED_PCA_ARTIFACTS = (
    'QC_Adjusted_Transformed_Counts.tsv',
    'PCA_Coordinates.tsv',
    'PCA_Explained_Variance.tsv',
    'PCA_Plot.png',
    'Adjusted_PCA_Coordinates.tsv',
    'Adjusted_PCA_Explained_Variance.tsv',
    'Adjusted_PCA_Plot.png',
)

def _get_transformed_counts(
    counts_T,
    metadata,
    design,
    continuous_factors=None,
    transform='vst',
    use_design=False,
    n_cpus=4,
):
    if transform == 'vst':
        dds = DeseqDataSet(
            counts=counts_T,
            metadata=metadata,
            design=design,
            continuous_factors=continuous_factors,
            n_cpus=n_cpus,
            quiet=True,
        )
        dds.vst(use_design=use_design)
        return pd.DataFrame(dds.layers["vst_counts"], index=counts_T.index, columns=counts_T.columns)

    if transform == 'log1p':
        return np.log1p(counts_T)

    raise QCValidationError(
        "The requested QC transform is not supported.",
        details={
            "reason": "unsupported_qc_transform",
            "transform": str(transform),
            "supported": ["vst", "log1p"],
        },
    )


def run_qc(
    counts_T,
    metadata,
    design_col,
    out_dir,
    design=None,
    continuous_factors=None,
    transform='vst',
    adjust_factors=None,
    use_design=False,
    label_samples=False,
    n_cpus=4,
    pca_top_n=500,
    biology_factors=None,
):
    """Generate experimental visualization-only QC artifacts.

    PCA genes are selected once from the unadjusted transformed matrix and the
    identical stable gene IDs are reused for ordinary and adjusted PCA.  These
    coordinates are not inputs to differential-expression testing.
    """
    nuisance_factors = (
        [adjust_factors]
        if isinstance(adjust_factors, str)
        else list(adjust_factors or [])
    )
    requested_biology_factors = (
        [biology_factors]
        if isinstance(biology_factors, str)
        else list(biology_factors or [design_col])
    )
    continuous_factor_names = (
        [continuous_factors]
        if isinstance(continuous_factors, str)
        else list(continuous_factors or [])
    )
    if nuisance_factors:
        protected_biology_factors = list(
            resolve_adjustment_biology_factors(
                design or f"~ {design_col}",
                design_col,
                requested_biology_factors,
                nuisance_factors,
            )
        )
    else:
        protected_biology_factors = [design_col]
        for factor in requested_biology_factors:
            if factor not in protected_biology_factors:
                protected_biology_factors.append(factor)
    os.makedirs(out_dir, exist_ok=True)
    metadata = _align_qc_metadata(counts_T.index, metadata)
    transformed = _get_transformed_counts(
        counts_T,
        metadata,
        design=design or f"~ {design_col}",
        continuous_factors=continuous_factors,
        transform=transform,
        use_design=use_design,
        n_cpus=n_cpus,
    )
    metadata = _align_qc_metadata(transformed.index, metadata)
    selection = select_top_variable_features(
        transformed.to_numpy(dtype=float),
        transformed.columns,
        pca_top_n,
    )
    if len(selection.names) < 2:
        raise QCValidationError(
            "QC PCA requires at least two positive-variance genes.",
            details={
                "reason": "insufficient_pca_dimensions",
                "selected_gene_count": len(selection.names),
            },
        )
    _guard_corrected_pca_outputs(out_dir)
    adjusted, adjustment_design = _residualize_covariates(
        transformed,
        metadata,
        nuisance_factors,
        protected_biology_factors,
        continuous_factor_names,
        return_design=True,
    )
    _write_pca_evidence(
        out_dir,
        selection,
        transform=transform,
        requested_top_n=pca_top_n,
        adjust_factors=nuisance_factors,
        biology_factors=protected_biology_factors,
        continuous_factors=continuous_factor_names,
        adjustment_design=adjustment_design,
    )
    transformed.to_csv(os.path.join(out_dir, 'QC_Transformed_Counts.tsv'), sep='\t')

    sample_qc = pd.DataFrame(index=counts_T.index)
    sample_qc['library_size'] = counts_T.sum(axis=1)
    sample_qc['detected_genes'] = (counts_T > 0).sum(axis=1)
    sample_qc = sample_qc.join(metadata)
    sample_qc.to_csv(os.path.join(out_dir, 'Sample_QC_Metrics.tsv'), sep='\t')

    plt.figure(figsize=(10, 5))
    sns.barplot(data=sample_qc.reset_index(), x=sample_qc.index.name or 'index', y='library_size', hue=design_col)
    plt.xticks(rotation=90)
    plt.title("Library Size by Sample")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Library_Size.png'), dpi=300)
    plt.close()

    plt.figure(figsize=(10, 5))
    sns.barplot(data=sample_qc.reset_index(), x=sample_qc.index.name or 'index', y='detected_genes', hue=design_col)
    plt.xticks(rotation=90)
    plt.title("Detected Genes by Sample")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Detected_Genes.png'), dpi=300)
    plt.close()
    
    selected_genes = list(selection.names)
    pca_df, _ = _save_pca_outputs(
        transformed.loc[:, selected_genes],
        metadata,
        design_col,
        out_dir,
        prefix='PCA',
        title_prefix=f'PCA ({transform.upper()})',
        label_samples=label_samples,
    )
    
    # Correlation
    corr = transformed.T.corr()
    corr.to_csv(os.path.join(out_dir, 'Sample_Correlation.tsv'), sep='\t')
    plt.figure(figsize=(10,8))
    sns.heatmap(corr, cmap='RdBu_r', center=0.8) # Red-Blue color map
    plt.title(f"Sample Correlation ({transform.upper()})")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Sample_Correlation.png'), dpi=300)
    plt.close()

    distances = pd.DataFrame(
        pairwise_distances(transformed, metric='euclidean'),
        index=transformed.index,
        columns=transformed.index,
    )
    distances.to_csv(os.path.join(out_dir, 'Sample_Distance.tsv'), sep='\t')
    plt.figure(figsize=(10, 8))
    sns.heatmap(distances, cmap='viridis')
    plt.title(f"Sample Distance ({transform.upper()})")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Sample_Distance.png'), dpi=300)
    plt.close()

    association_rows = []
    pcs = pca_df[['PC1', 'PC2']]
    for col in metadata.columns:
        series = metadata[col]
        if col in continuous_factor_names:
            numeric_series = pd.to_numeric(series, errors='raise')
            assoc_pc1 = abs(numeric_series.corr(pcs['PC1']))
            assoc_pc2 = abs(numeric_series.corr(pcs['PC2']))
            association_rows.append({
                'factor': col,
                'type': 'continuous',
                'score_pc1': float(assoc_pc1) if pd.notna(assoc_pc1) else np.nan,
                'score_pc2': float(assoc_pc2) if pd.notna(assoc_pc2) else np.nan,
                'best_score': np.nanmax([assoc_pc1, assoc_pc2]),
            })
        else:
            labels = series.astype(str)
            if labels.nunique() < 2 or labels.nunique() >= len(labels):
                continue
            score = silhouette_score(pcs, labels)
            association_rows.append({
                'factor': col,
                'type': 'categorical',
                'score_pc1': np.nan,
                'score_pc2': np.nan,
                'best_score': float(score),
            })
    if association_rows:
        assoc_df = pd.DataFrame(association_rows).sort_values('best_score', ascending=False)
        assoc_df.to_csv(os.path.join(out_dir, 'QC_Metadata_Associations.tsv'), sep='\t', index=False)
    else:
        pd.DataFrame(columns=['factor', 'type', 'score_pc1', 'score_pc2', 'best_score']).to_csv(
            os.path.join(out_dir, 'QC_Metadata_Associations.tsv'),
            sep='\t',
            index=False,
        )

    adjusted.to_csv(os.path.join(out_dir, 'QC_Adjusted_Transformed_Counts.tsv'), sep='\t')
    _save_pca_outputs(
        adjusted.loc[:, selected_genes],
        metadata,
        design_col,
        out_dir,
        prefix='Adjusted_PCA',
        title_prefix='Adjusted PCA',
        label_samples=label_samples,
    )


def _align_qc_metadata(sample_ids, metadata):
    reorder = metadata_reorder_indices(sample_ids, metadata.index)
    aligned = metadata.iloc[reorder].copy()
    aligned.index = pd.Index(sample_ids, name=metadata.index.name)
    return aligned


def _guard_corrected_pca_outputs(out_dir):
    existing = [
        name
        for name in _CORRECTED_PCA_ARTIFACTS
        if os.path.exists(os.path.join(out_dir, name))
    ]
    method_path = os.path.join(out_dir, _PCA_METHOD_FILE)
    selected_path = os.path.join(out_dir, _PCA_SELECTED_GENES_FILE)
    if not existing and not os.path.exists(method_path):
        return
    if not os.path.exists(method_path):
        raise QCValidationError(
            "Refusing to overwrite legacy PCA artifacts without method provenance.",
            details={
                "reason": "legacy_pca_output_conflict",
                "existing_artifacts": existing,
                "required_provenance": method_path,
                "recovery": "Archive or remove the old QC directory before rerunning.",
            },
        )
    try:
        with open(method_path, 'r', encoding='utf-8') as handle:
            provenance = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise QCValidationError(
            "Refusing to overwrite PCA artifacts because provenance is unreadable.",
            details={
                "reason": "pca_provenance_unreadable",
                "provenance_path": method_path,
                "cause_type": type(exc).__name__,
                "cause_message": str(exc),
            },
            cause=exc,
        ) from exc
    if provenance.get('method_id') != PCA_METHOD_ID:
        raise QCValidationError(
            "Refusing to overwrite PCA artifacts from a different or legacy method.",
            details={
                "reason": "pca_method_conflict",
                "provenance_path": method_path,
                "observed_method_id": provenance.get('method_id'),
                "required_method_id": PCA_METHOD_ID,
            },
        )
    if existing and not os.path.exists(selected_path):
        raise QCValidationError(
            "Refusing to overwrite PCA artifacts without selected-gene evidence.",
            details={
                "reason": "pca_selection_evidence_missing",
                "required_evidence": selected_path,
                "existing_artifacts": existing,
            },
        )


def _write_pca_evidence(
    out_dir,
    selection,
    *,
    transform,
    requested_top_n,
    adjust_factors,
    biology_factors,
    continuous_factors,
    adjustment_design,
):
    selected_rows = pd.DataFrame({
        'rank': np.arange(1, len(selection.names) + 1, dtype=int),
        'gene_id': selection.names,
        'sample_variance_unadjusted': selection.variances,
    })
    selected_rows.to_csv(
        os.path.join(out_dir, _PCA_SELECTED_GENES_FILE),
        sep='\t',
        index=False,
    )
    selected_digest = hashlib.sha256(
        ''.join(f'{gene_id}\n' for gene_id in selection.names).encode('utf-8')
    ).hexdigest()
    continuous = set(continuous_factors)
    factor_semantics = [
        {
            'factor': factor,
            'role': role,
            'kind': 'numeric' if factor in continuous else 'categorical',
        }
        for role, factors in (
            ('biology', biology_factors),
            ('nuisance', adjust_factors),
        )
        for factor in factors
    ]
    provenance = {
        'schema_version': '1.0',
        'status': 'experimental',
        'analysis_role': 'visualization_only',
        'method_id': PCA_METHOD_ID,
        'transform': str(transform),
        'selection_source': 'unadjusted_transformed_expression',
        'selection_rule': 'top_positive_sample_variance',
        'selection_variance_ddof': 1,
        'requested_top_n': int(requested_top_n),
        'selected_gene_count': len(selection.names),
        'selected_gene_sha256': selected_digest,
        'ordinary_and_adjusted_pca_reuse_exact_gene_ids': True,
        'pca_center': True,
        'pca_scale': False,
        'pca_components': 2,
        'pca_implementation': 'numpy.linalg.svd',
        'biology_factors': list(biology_factors),
        'nuisance_factors': list(adjust_factors),
        'factor_semantics': factor_semantics,
        'adjustment_applied': bool(adjust_factors),
        'biology_design_columns': (
            list(adjustment_design.biology_columns)
            if adjustment_design is not None
            else []
        ),
        'nuisance_design_columns': (
            list(adjustment_design.nuisance_columns)
            if adjustment_design is not None
            else []
        ),
        'joint_design_rank': (
            int(adjustment_design.rank) if adjustment_design is not None else None
        ),
        'joint_design_residual_df': (
            int(adjustment_design.residual_df)
            if adjustment_design is not None
            else None
        ),
        'categorical_biology_encoding': 'treatment_contrasts',
        'categorical_nuisance_encoding': 'limma_style_sum_contrasts',
        'numeric_nuisance_encoding': 'sample_mean_centered',
        'adjustment_model': (
            'Y = X_biology B + X_nuisance Gamma + E; '
            'subtract X_nuisance Gamma only'
        ),
        'adjusted_values_enter_differential_expression': False,
    }
    with open(os.path.join(out_dir, _PCA_METHOD_FILE), 'w', encoding='utf-8') as handle:
        json.dump(provenance, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def _residualize_covariates(
    transformed,
    metadata,
    adjust_factors,
    biology_factors,
    continuous_factors,
    *,
    return_design=False,
):
    """Remove nuisance terms from a joint model for visualization only."""

    metadata = _align_qc_metadata(transformed.index, metadata)
    factors = list(adjust_factors)
    if not factors:
        unchanged = transformed.copy()
        return (unchanged, None) if return_design else unchanged
    biology = (
        [biology_factors]
        if isinstance(biology_factors, str)
        else list(biology_factors)
    )
    continuous = set(continuous_factors)
    validate_factor_names(metadata.columns, biology, factors)

    def factor_spec(factor):
        series = metadata[factor]
        if series.isna().any():
            raise QCValidationError(
                f"Factor '{factor}' contains missing values.",
                details={"reason": "missing_factor_value", "factor": factor},
            )
        is_numeric = factor in continuous
        kind = 'numeric' if is_numeric else 'categorical'
        values = series.astype(object).to_numpy()
        return FactorSpec(factor, values, kind)

    design = build_joint_design(
        [factor_spec(factor) for factor in biology],
        [factor_spec(factor) for factor in factors],
    )
    adjusted = remove_nuisance_effects(
        transformed.to_numpy(dtype=float),
        design,
    )
    adjusted_frame = pd.DataFrame(
        adjusted,
        index=transformed.index.copy(),
        columns=transformed.columns.copy(),
    )
    return (adjusted_frame, design) if return_design else adjusted_frame


def _save_pca_outputs(
    transformed,
    metadata,
    design_col,
    out_dir,
    prefix,
    title_prefix,
    label_samples=False,
):
    metadata = _align_qc_metadata(transformed.index, metadata)
    pca = centered_unscaled_pca(transformed.to_numpy(dtype=float), n_components=2)
    pca_df = pd.DataFrame(pca.coordinates, index=metadata.index, columns=['PC1', 'PC2'])
    pca_df = pd.concat([pca_df, metadata], axis=1)
    pca_df.to_csv(os.path.join(out_dir, f'{prefix}_Coordinates.tsv'), sep='\t')
    pd.DataFrame({
        'component': ['PC1', 'PC2'],
        'explained_variance_ratio': pca.explained_variance_ratio,
    }).to_csv(os.path.join(out_dir, f'{prefix}_Explained_Variance.tsv'), sep='\t', index=False)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=design_col, style=design_col, s=150)
    if label_samples:
        for sample_id, row in pca_df[['PC1', 'PC2']].iterrows():
            plt.text(row['PC1'], row['PC2'], str(sample_id), fontsize=8)
    plt.title(
        f"{title_prefix}: "
        f"PC1 {pca.explained_variance_ratio[0]:.1%}, "
        f"PC2 {pca.explained_variance_ratio[1]:.1%}"
    )
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{prefix}_Plot.png'), dpi=300)
    plt.close()
    return pca_df, pca

def fit_model(counts_T, metadata, design, continuous_factors=None, n_cpus=4):
    print("  Fitting DESeq2 model...")
    dds = DeseqDataSet(
        counts=counts_T, 
        metadata=metadata, 
        design=design,
        continuous_factors=continuous_factors,
        n_cpus=n_cpus, quiet=True
    )
    dds.deseq2()
    return dds

def _shrinkable_coeff_name(stat_res, contrast):
    candidate = f"{contrast['factor']}[T.{contrast['treatment']}]"
    return candidate if candidate in stat_res.LFC.columns else None


def run_contrasts(
    dds,
    contrasts,
    default_factor,
    out_dir,
    n_cpus=4,
    reference_levels=None,
    shrink_lfc=True,
    shrink_adapt=True,
):
    os.makedirs(out_dir, exist_ok=True)
    results = {}
    
    for contrast in contrasts:
        name = data.contrast_name(contrast, default_factor)
        factor = contrast['factor']
        treat = contrast['treatment']
        ctrl = contrast['control']
        print(f"  Comparing: {name}")
        try:
            stat_res = DeseqStats(dds, contrast=[factor, treat, ctrl], n_cpus=n_cpus, quiet=True)
            stat_res.summary()
            df = stat_res.results_df.copy()
            df['contrast_factor'] = factor
            df['contrast_treatment'] = treat
            df['contrast_control'] = ctrl
            df['lfc_shrink_applied'] = False

            ref_for_factor = (reference_levels or {}).get(factor)
            shrink_coeff = _shrinkable_coeff_name(stat_res, contrast)
            can_shrink = shrink_lfc and ctrl == ref_for_factor and shrink_coeff is not None

            if can_shrink:
                df['log2FoldChange_raw'] = df['log2FoldChange']
                df['lfcSE_raw'] = df['lfcSE']
                stat_res.lfc_shrink(coeff=shrink_coeff, adapt=shrink_adapt)
                stat_res.summary()
                shrunk = stat_res.results_df
                df['log2FoldChange'] = shrunk['log2FoldChange']
                df['lfcSE'] = shrunk['lfcSE']
                df['log2FoldChange_shrunk'] = shrunk['log2FoldChange']
                df['lfcSE_shrunk'] = shrunk['lfcSE']
                df['lfc_shrink_applied'] = True

            df.to_csv(os.path.join(out_dir, f"{name}.csv"))
            results[name] = df
        except Exception as e:
            print(f"    [Warning] Failed contrast {name}: {e}")
            
    return results


def write_contrast_summary(results, out_dir, padj_thresh=0.05, logfc_thresh=1.0):
    summary_rows = []
    for name, df in results.items():
        clean = df.dropna(subset=['padj', 'log2FoldChange'])
        summary_rows.append({
            'contrast': name,
            'n_tested_genes': int(clean.shape[0]),
            'n_sig': int(((clean['padj'] < padj_thresh)).sum()),
            'n_up': int(((clean['padj'] < padj_thresh) & (clean['log2FoldChange'] > logfc_thresh)).sum()),
            'n_down': int(((clean['padj'] < padj_thresh) & (clean['log2FoldChange'] < -logfc_thresh)).sum()),
        })

    summary_df = pd.DataFrame(summary_rows).sort_values('contrast')
    out_path = os.path.join(out_dir, "_contrast_summary.csv")
    summary_df.to_csv(out_path, index=False)
    return out_path

def _top_label_genes(df, n=10):
    clean = df.dropna(subset=['padj', 'log2FoldChange']).copy()
    if clean.empty or n <= 0:
        return clean.iloc[0:0]
    clean['rank_score'] = (-np.log10(clean['padj'].clip(lower=1e-300))) * clean['log2FoldChange'].abs()
    return clean.sort_values(['rank_score', 'padj'], ascending=[False, True]).head(n)


def plot_volcano(results, out_dir, padj_thresh=0.05, logfc_thresh=1.0, top_label_genes=10):
    os.makedirs(out_dir, exist_ok=True)
    for name, df in results.items():
        df = df.dropna(subset=['padj', 'log2FoldChange']).copy()
        df['padj_plot'] = df['padj'].clip(lower=1e-300)
        
        # Vectorized color assignment
        # Red: Sig Up
        mask_up = (df['padj'] < padj_thresh) & (df['log2FoldChange'] > logfc_thresh)
        # Blue: Sig Down
        mask_down = (df['padj'] < padj_thresh) & (df['log2FoldChange'] < -logfc_thresh)
        
        plt.figure(figsize=(6,5))
        plt.scatter(df['log2FoldChange'], -np.log10(df['padj_plot']), c='lightgrey', s=5, alpha=0.5)
        plt.scatter(df.loc[mask_up, 'log2FoldChange'], -np.log10(df.loc[mask_up, 'padj_plot']), c='red', s=10, label='Up')
        plt.scatter(df.loc[mask_down, 'log2FoldChange'], -np.log10(df.loc[mask_down, 'padj_plot']), c='blue', s=10, label='Down')
        
        plt.axhline(-np.log10(padj_thresh), ls='--', c='k', lw=0.5)
        plt.axvline(logfc_thresh, ls='--', c='k', lw=0.5)
        plt.axvline(-logfc_thresh, ls='--', c='k', lw=0.5)
        
        plt.title(name)
        plt.xlabel("log2 Fold Change")
        plt.ylabel("-log10(padj)")
        for gene, row in _top_label_genes(df, top_label_genes).iterrows():
            plt.text(row['log2FoldChange'], -np.log10(max(row['padj_plot'], 1e-300)), str(gene), fontsize=7)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"Volcano_{name}.png"), dpi=300)
        plt.close()


def plot_ma(results, out_dir, padj_thresh=0.05, logfc_thresh=1.0, top_label_genes=10):
    os.makedirs(out_dir, exist_ok=True)
    for name, df in results.items():
        clean = df.dropna(subset=['baseMean', 'log2FoldChange', 'padj']).copy()
        sig = clean['padj'] < padj_thresh
        strong = sig & (clean['log2FoldChange'].abs() > logfc_thresh)

        plt.figure(figsize=(6, 5))
        plt.scatter(np.log10(clean['baseMean'] + 1), clean['log2FoldChange'], c='lightgrey', s=6, alpha=0.5)
        plt.scatter(
            np.log10(clean.loc[strong, 'baseMean'] + 1),
            clean.loc[strong, 'log2FoldChange'],
            c='crimson',
            s=8,
            alpha=0.7,
        )
        plt.axhline(0, ls='--', c='k', lw=0.5)
        plt.axhline(logfc_thresh, ls='--', c='k', lw=0.5)
        plt.axhline(-logfc_thresh, ls='--', c='k', lw=0.5)
        plt.xlabel('log10(baseMean + 1)')
        plt.ylabel('log2 Fold Change')
        plt.title(f"MA Plot: {name}")
        for gene, row in _top_label_genes(clean, top_label_genes).iterrows():
            plt.text(np.log10(row['baseMean'] + 1), row['log2FoldChange'], str(gene), fontsize=7)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"MA_{name}.png"), dpi=300)
        plt.close()


def write_cooks_diagnostics(dds, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    cooks_matrix = pd.DataFrame(dds.layers["cooks"], index=dds.obs_names, columns=dds.var_names)
    cooks_matrix.to_csv(os.path.join(out_dir, "_cooks_matrix.tsv"), sep='\t')

    gene_report = pd.DataFrame(index=dds.var_names)
    gene_report['max_cooks'] = cooks_matrix.max(axis=0)
    gene_report['mean_cooks'] = cooks_matrix.mean(axis=0)
    gene_report['cooks_outlier'] = pd.Series(dds.cooks_outlier(), index=dds.var_names).astype(bool)
    if 'replaced' in dds.var:
        gene_report['refit_replaced'] = dds.var['replaced'].astype(bool).values
    gene_report.to_csv(os.path.join(out_dir, "_cooks_gene_report.tsv"), sep='\t')

    sample_report = pd.DataFrame(index=dds.obs_names)
    sample_report['max_cooks'] = cooks_matrix.max(axis=1)
    sample_report['mean_cooks'] = cooks_matrix.mean(axis=1)
    sample_report['n_genes_above_99pct'] = (cooks_matrix.gt(cooks_matrix.quantile(0.99).fillna(np.inf), axis=1)).sum(axis=1)
    sample_report.to_csv(os.path.join(out_dir, "_cooks_sample_report.tsv"), sep='\t')

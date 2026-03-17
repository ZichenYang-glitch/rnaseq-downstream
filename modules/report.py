import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from modules import data


def create_master_table(
    metadata,
    design_col,
    results_dict,
    tpm_file,
    counts_T,
    out_dir,
    strip_gene_version=True,
    annotation_df=None,
):
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Base Expression (TPM or Counts)
    if os.path.exists(tpm_file):
        base = data.load_expression_matrix(
            tpm_file,
            metadata.index,
            strip_gene_version=strip_gene_version,
            aggregate='mean',
        )
    else:
        base = counts_T.T.copy()
        
    # Calculate group means
    groups = metadata[design_col].unique()
    for g in groups:
        samples = metadata[metadata[design_col] == g].index
        valid = [s for s in samples if s in base.columns]
        if valid:
            base[f"Mean_{g}"] = base[valid].mean(axis=1)
            
    # 2. Add Stats
    final = base.copy()
    final.index.name = 'gene_id'
    if annotation_df is not None:
        final = final.join(annotation_df, how='left')
    for name, df in results_dict.items():
        stat_cols = [col for col in [
            'log2FoldChange',
            'padj',
            'pvalue',
            'stat',
            'log2FoldChange_raw',
            'log2FoldChange_shrunk',
            'lfc_shrink_applied',
        ] if col in df.columns]
        sub = df[stat_cols].copy()
        sub.columns = [f"{col}_{name}" for col in stat_cols]
        final = final.join(sub, how='left')
        
    out_path = os.path.join(out_dir, "Master_Expression_Table.csv")
    final.to_csv(out_path)
    print(f"  Master table saved to: {out_path}")


def create_analysis_summary(
    metadata,
    design_col,
    results_dict,
    qc_dir,
    validation_dir,
    out_dir,
    upstream_provenance=None,
):
    os.makedirs(out_dir, exist_ok=True)

    lines = [
        "# Analysis Summary",
        "",
    ]

    if upstream_provenance:
        lines.extend([
            "## Upstream Provenance",
            "",
        ])
        for key, value in upstream_provenance.items():
            lines.append(f"- {key}: {value}")
        lines.extend(["",])

    lines.extend([
        "## Samples",
        "",
        f"- Total samples: {metadata.shape[0]}",
        f"- Design factor: `{design_col}`",
        "",
        "### Group Counts",
        "",
    ])

    for group, n in metadata[design_col].astype(str).value_counts().sort_index().items():
        lines.append(f"- {group}: {n}")

    lines.extend(["", "## QC", ""])
    pca_var_path = os.path.join(qc_dir, "PCA_Explained_Variance.tsv")
    assoc_path = os.path.join(qc_dir, "QC_Metadata_Associations.tsv")
    sample_qc_path = os.path.join(qc_dir, "Sample_QC_Metrics.tsv")
    if os.path.exists(sample_qc_path):
        sample_qc = pd.read_csv(sample_qc_path, sep='\t', index_col=0)
        lines.append(
            f"- Library size range: {int(sample_qc['library_size'].min()):,} to "
            f"{int(sample_qc['library_size'].max()):,}"
        )
        lines.append(
            f"- Detected genes range: {int(sample_qc['detected_genes'].min()):,} to "
            f"{int(sample_qc['detected_genes'].max()):,}"
        )
    if os.path.exists(pca_var_path):
        pca_var = pd.read_csv(pca_var_path, sep='\t')
        if pca_var.shape[0] >= 2:
            lines.append(
                f"- PCA explained variance: PC1 {pca_var.loc[0, 'explained_variance_ratio']:.2%}, "
                f"PC2 {pca_var.loc[1, 'explained_variance_ratio']:.2%}"
            )
    if os.path.exists(assoc_path):
        assoc = pd.read_csv(assoc_path, sep='\t')
        if not assoc.empty:
            top = assoc.sort_values('best_score', ascending=False).iloc[0]
            lines.append(
                f"- Strongest QC-associated metadata field: `{top['factor']}` "
                f"({top['type']}, score={top['best_score']:.3f})"
            )

    lines.extend(["", "## Differential Expression", ""])
    contrast_rows = []
    for name, df in results_dict.items():
        clean = df.dropna(subset=['padj', 'log2FoldChange'])
        contrast_rows.append({
            'contrast': name,
            'tested': int(clean.shape[0]),
            'significant': int((clean['padj'] < 0.05).sum()),
            'up': int(((clean['padj'] < 0.05) & (clean['log2FoldChange'] > 1)).sum()),
            'down': int(((clean['padj'] < 0.05) & (clean['log2FoldChange'] < -1)).sum()),
            'shrinkage': bool(df.get('lfc_shrink_applied', pd.Series([False])).any()),
        })
    if contrast_rows:
        summary_df = pd.DataFrame(contrast_rows).sort_values('significant', ascending=False)
        for _, row in summary_df.iterrows():
            lines.append(
                f"- {row['contrast']}: tested={row['tested']}, sig={row['significant']}, "
                f"up={row['up']}, down={row['down']}, shrinkage={row['shrinkage']}"
            )

    validation_file = os.path.join(validation_dir, "validated_inputs.txt")
    if os.path.exists(validation_file):
        lines.extend(["", "## Validation", ""])
        with open(validation_file, 'r', encoding='utf-8') as handle:
            for line in handle:
                key, value = line.rstrip('\n').split('\t', 1)
                lines.append(f"- {key}: {value}")

    out_path = os.path.join(out_dir, "Analysis_Summary.md")
    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write("\n".join(lines) + "\n")
    print(f"  Analysis summary saved to: {out_path}")


def _load_transformed_counts(qc_dir):
    path = os.path.join(qc_dir, "QC_Transformed_Counts.tsv")
    return pd.read_csv(path, sep='\t', index_col=0)


def _zscore_rows(df):
    centered = df.sub(df.mean(axis=1), axis=0)
    scaled = centered.div(df.std(axis=1).replace(0, np.nan), axis=0)
    return scaled.fillna(0)


def _save_placeholder_plot(path, title, message):
    plt.figure(figsize=(8, 4))
    plt.axis('off')
    plt.title(title)
    plt.text(0.5, 0.5, message, ha='center', va='center')
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def create_heatmap_summaries(
    metadata,
    design_col,
    results_dict,
    qc_dir,
    out_dir,
    top_variable_genes=50,
    top_de_genes_per_contrast=30,
    annotation_factors=None,
):
    os.makedirs(out_dir, exist_ok=True)
    transformed = _load_transformed_counts(qc_dir)
    ordered_samples = metadata.sort_values(design_col).index.tolist()
    transformed = transformed.loc[ordered_samples]

    variable_genes = transformed.var(axis=0).sort_values(ascending=False).head(top_variable_genes).index.tolist()
    variable_path = os.path.join(out_dir, "Top_Variable_Genes_Heatmap.png")
    if variable_genes:
        var_mat = _zscore_rows(transformed[variable_genes].T)
        sns.clustermap(var_mat, cmap='vlag', col_cluster=False, figsize=(10, 10))
        plt.savefig(variable_path, dpi=300)
        plt.close('all')
    else:
        _save_placeholder_plot(variable_path, "Top Variable Genes Heatmap", "No genes available.")

    de_genes = []
    for _, df in results_dict.items():
        clean = df.dropna(subset=['padj', 'log2FoldChange']).copy()
        ranked = clean.sort_values(['padj', 'log2FoldChange'], ascending=[True, False])
        de_genes.extend(ranked.head(top_de_genes_per_contrast).index.tolist())
    de_genes = list(dict.fromkeys([gene for gene in de_genes if gene in transformed.columns]))
    de_path = os.path.join(out_dir, "Top_DE_Genes_Heatmap.png")
    if de_genes:
        de_mat = _zscore_rows(transformed[de_genes].T)
        sns.clustermap(de_mat, cmap='RdBu_r', col_cluster=False, figsize=(10, 12))
        plt.savefig(de_path, dpi=300)
        plt.close('all')
    else:
        _save_placeholder_plot(de_path, "Top DE Genes Heatmap", "No DE genes available.")

    cluster_path = os.path.join(out_dir, "Annotated_Top_Variable_Clustermap.png")
    annotation_factors = [factor for factor in (annotation_factors or [design_col]) if factor in metadata.columns]
    if variable_genes:
        annot = metadata.loc[ordered_samples, annotation_factors].astype(str)
        palettes = {}
        color_frames = []
        for col in annotation_factors:
            levels = sorted(annot[col].unique())
            palette = dict(zip(levels, sns.color_palette('tab10', n_colors=max(3, len(levels)))[:len(levels)]))
            palettes[col] = palette
            color_frames.append(annot[col].map(palette).rename(col))
        col_colors = pd.concat(color_frames, axis=1) if color_frames else None
        g = sns.clustermap(var_mat, cmap='vlag', col_cluster=True, row_cluster=True, col_colors=col_colors, figsize=(12, 10))
        g.fig.suptitle("Annotated Top Variable Genes Clustermap", y=1.02)
        g.savefig(cluster_path, dpi=300)
        plt.close('all')
    else:
        _save_placeholder_plot(cluster_path, "Annotated Clustermap", "No genes available.")


def create_deg_summary_plot(results_dict, out_dir, padj_thresh=0.05, logfc_thresh=1.0):
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    for name, df in results_dict.items():
        clean = df.dropna(subset=['padj', 'log2FoldChange'])
        rows.append({
            'contrast': name,
            'up': int(((clean['padj'] < padj_thresh) & (clean['log2FoldChange'] > logfc_thresh)).sum()),
            'down': int(((clean['padj'] < padj_thresh) & (clean['log2FoldChange'] < -logfc_thresh)).sum()),
        })
    summary = pd.DataFrame(rows).sort_values('contrast')
    summary.to_csv(os.path.join(out_dir, "DEG_Counts.tsv"), sep='\t', index=False)

    plt.figure(figsize=(10, 6))
    plt.bar(summary['contrast'], summary['up'], color='crimson', label='Up')
    plt.bar(summary['contrast'], -summary['down'], color='steelblue', label='Down')
    plt.axhline(0, color='black', lw=0.8)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Gene Count")
    plt.title("Differentially Expressed Genes by Contrast")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "DEG_Counts_Barplot.png"), dpi=300)
    plt.close()


def create_sample_outlier_report(qc_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    metrics = pd.read_csv(os.path.join(qc_dir, "Sample_QC_Metrics.tsv"), sep='\t', index_col=0)
    distances = pd.read_csv(os.path.join(qc_dir, "Sample_Distance.tsv"), sep='\t', index_col=0)

    mean_distance = distances.replace(0, np.nan).mean(axis=1)
    z_mean_distance = (mean_distance - mean_distance.mean()) / (mean_distance.std(ddof=0) or 1.0)
    z_library_size = (metrics['library_size'] - metrics['library_size'].mean()) / (metrics['library_size'].std(ddof=0) or 1.0)
    z_detected = (metrics['detected_genes'] - metrics['detected_genes'].mean()) / (metrics['detected_genes'].std(ddof=0) or 1.0)

    report = metrics.copy()
    report['mean_sample_distance'] = mean_distance
    report['z_mean_sample_distance'] = z_mean_distance
    report['z_library_size'] = z_library_size
    report['z_detected_genes'] = z_detected
    report['potential_outlier'] = (
        report['z_mean_sample_distance'].abs().gt(2)
        | report['z_library_size'].abs().gt(2)
        | report['z_detected_genes'].abs().gt(2)
    )
    report.to_csv(os.path.join(out_dir, "Sample_Outlier_Report.tsv"), sep='\t')


def _collect_gsea_tables(gsea_dir):
    rows = []
    for path in glob.glob(os.path.join(gsea_dir, "*", "*", "*.csv")):
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if 'Term' not in df.columns or 'NES' not in df.columns:
            continue
        fdr_col = next((col for col in ['FDR q-val', 'FDR', 'fdr'] if col in df.columns), None)
        if fdr_col is None:
            continue
        parts = path.split(os.sep)
        contrast = parts[-3]
        gene_set = parts[-2]
        sub = df[['Term', 'NES', fdr_col]].copy()
        sub.columns = ['term', 'NES', 'FDR']
        sub['contrast'] = contrast
        sub['gene_set'] = gene_set
        rows.append(sub)
    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(columns=['term', 'NES', 'FDR', 'contrast', 'gene_set'])


def create_gsea_summary(gsea_dir, out_dir, fdr_thresh=0.25, top_terms=20):
    os.makedirs(out_dir, exist_ok=True)
    summary = _collect_gsea_tables(gsea_dir)
    summary.to_csv(os.path.join(out_dir, "GSEA_Summary.tsv"), sep='\t', index=False)

    plot_path = os.path.join(out_dir, "GSEA_Summary_Dotplot.png")
    if summary.empty:
        _save_placeholder_plot(plot_path, "GSEA Summary", "No GSEA results found.")
        return

    filtered = summary[summary['FDR'] <= fdr_thresh].copy()
    if filtered.empty:
        _save_placeholder_plot(plot_path, "GSEA Summary", "No pathways passed the FDR threshold.")
        return

    filtered['abs_nes'] = filtered['NES'].abs()
    top = filtered.sort_values(['abs_nes', 'FDR'], ascending=[False, True]).head(top_terms)
    top['neglog10_fdr'] = -np.log10(top['FDR'].clip(lower=1e-300))

    plt.figure(figsize=(12, max(4, 0.4 * len(top))))
    sns.scatterplot(
        data=top,
        x='contrast',
        y='term',
        size='neglog10_fdr',
        hue='NES',
        palette='coolwarm',
        sizes=(30, 250),
        edgecolor='black',
    )
    plt.title("Top Enriched Pathways Across Contrasts")
    plt.xlabel("Contrast")
    plt.ylabel("Pathway")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()

    matrix_path = os.path.join(out_dir, "GSEA_Term_Contrast_Matrix.tsv")
    heatmap_path = os.path.join(out_dir, "GSEA_Term_Clustermap.png")
    pivot = top.pivot_table(index='term', columns='contrast', values='NES', aggfunc='mean').fillna(0)
    pivot.to_csv(matrix_path, sep='\t')
    if not pivot.empty:
        g = sns.clustermap(pivot, cmap='coolwarm', center=0, figsize=(10, max(6, 0.3 * pivot.shape[0])))
        g.fig.suptitle("GSEA Term Clustering Across Contrasts", y=1.02)
        g.savefig(heatmap_path, dpi=300)
        plt.close('all')
    else:
        _save_placeholder_plot(heatmap_path, "GSEA Term Clustering", "No terms available.")


def create_report_index(output_dir, out_dir, upstream_provenance=None):
    os.makedirs(out_dir, exist_ok=True)
    lines = [
        "# Report Index",
        "",
        "## Core Outputs",
        "",
        "- `00_Validation/validated_inputs.txt`",
        "- `00_Validation/upstream_provenance.tsv`",
        "- `01_QC/PCA_Plot.png`",
        "- `01_QC/Adjusted_PCA_Plot.png`",
        "- `01_QC/QC_Metadata_Associations.tsv`",
        "- `02_DESeq2_Stats/_contrast_summary.csv`",
        "- `03_Volcano_Plots/` for volcano and MA plots",
        "- `05_Summary/Analysis_Summary.md`",
        "- `05_Summary/Annotated_Top_Variable_Clustermap.png`",
        "- `05_Summary/GSEA_Summary.tsv`",
        "- `05_Summary/Sample_Outlier_Report.tsv`",
        "",
        "## Navigation",
        "",
        f"- Output root: `{output_dir}`",
    ]
    if upstream_provenance:
        lines.extend(["", "## Upstream Provenance", ""])
        for key, value in upstream_provenance.items():
            lines.append(f"- {key}: {value}")
    out_path = os.path.join(out_dir, "Report_Index.md")
    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write("\n".join(lines) + "\n")


def create_qc_adjustment_comparison(qc_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    orig_var = pd.read_csv(os.path.join(qc_dir, "PCA_Explained_Variance.tsv"), sep='\t')
    adj_var = pd.read_csv(os.path.join(qc_dir, "Adjusted_PCA_Explained_Variance.tsv"), sep='\t')
    assoc = pd.read_csv(os.path.join(qc_dir, "QC_Metadata_Associations.tsv"), sep='\t')

    rows = [
        {
            'metric': 'PC1_explained_variance_before',
            'value': orig_var.loc[orig_var['component'] == 'PC1', 'explained_variance_ratio'].iloc[0],
        },
        {
            'metric': 'PC2_explained_variance_before',
            'value': orig_var.loc[orig_var['component'] == 'PC2', 'explained_variance_ratio'].iloc[0],
        },
        {
            'metric': 'PC1_explained_variance_after',
            'value': adj_var.loc[adj_var['component'] == 'PC1', 'explained_variance_ratio'].iloc[0],
        },
        {
            'metric': 'PC2_explained_variance_after',
            'value': adj_var.loc[adj_var['component'] == 'PC2', 'explained_variance_ratio'].iloc[0],
        },
    ]
    if not assoc.empty:
        top = assoc.sort_values('best_score', ascending=False).iloc[0]
        rows.append({
            'metric': 'top_metadata_association_before',
            'value': f"{top['factor']} ({top['best_score']:.3f})",
        })

    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "QC_Adjustment_Comparison.tsv"), sep='\t', index=False)


def create_html_report_index(output_dir, out_dir, upstream_provenance=None):
    os.makedirs(out_dir, exist_ok=True)
    contrast_summary_path = os.path.join(output_dir, "02_DESeq2_Stats", "_contrast_summary.csv")
    outlier_path = os.path.join(out_dir, "Sample_Outlier_Report.tsv")
    gsea_summary_path = os.path.join(out_dir, "GSEA_Summary.tsv")
    upstream_path = os.path.join(output_dir, "00_Validation", "upstream_provenance.tsv")

    def rel(path):
        return os.path.relpath(path, out_dir).replace(os.sep, "/")

    def preview_table(path, sep="\t", n=8):
        if not os.path.exists(path):
            return "<p>Not available.</p>"
        try:
            df = pd.read_csv(path, sep=sep, engine='python').head(n)
        except Exception:
            return f'<p><a href="{rel(path)}">Open file</a></p>'
        if df.empty:
            return "<p>No rows available.</p>"
        headers = "".join(f"<th>{col}</th>" for col in df.columns)
        rows = []
        for _, row in df.iterrows():
            rows.append("<tr>" + "".join(f"<td>{row[col]}</td>" for col in df.columns) + "</tr>")
        return (
            '<div class="table-wrap"><table><thead><tr>'
            + headers
            + "</tr></thead><tbody>"
            + "".join(rows)
            + "</tbody></table></div>"
        )

    volcano_dir = os.path.join(output_dir, "03_Volcano_Plots")
    contrast_cards = []
    for volcano_path in sorted(glob.glob(os.path.join(volcano_dir, "Volcano_*.png"))):
        name = os.path.basename(volcano_path).replace("Volcano_", "").replace(".png", "")
        ma_path = os.path.join(volcano_dir, f"MA_{name}.png")
        stat_path = os.path.join(output_dir, "02_DESeq2_Stats", f"{name}.csv")
        detail_path = os.path.join(out_dir, f"contrast_{name}.html")
        contrast_cards.append(
            f"""
            <div class="card contrast-card">
              <h3>{name}</h3>
              <p><a href="{rel(detail_path)}">Detail page</a> | <a href="{rel(stat_path)}">DE table</a> | <a href="{rel(volcano_path)}">Volcano</a> | <a href="{rel(ma_path)}">MA</a></p>
              <a href="{rel(volcano_path)}"><img class="thumb" src="{rel(volcano_path)}" alt="Volcano {name}"></a>
            </div>
            """
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>RNA-seq Analysis Portal</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; line-height: 1.5; color: #222; background: #fafafa; }}
    h1, h2, h3 {{ margin-bottom: 0.5rem; }}
    ul {{ margin-top: 0; padding-left: 1.2rem; }}
    code {{ background: #f3f3f3; padding: 0.1rem 0.3rem; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 1rem; }}
    .card {{ border: 1px solid #ddd; border-radius: 10px; padding: 1rem; background: #fff; box-shadow: 0 1px 3px rgba(0,0,0,0.05); }}
    .thumb {{ width: 100%; border: 1px solid #eee; border-radius: 6px; margin-top: 0.75rem; }}
    a {{ color: #0b57d0; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .hero {{ margin-bottom: 1.5rem; }}
    .section {{ margin-top: 2rem; }}
    .table-wrap {{ overflow-x: auto; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.92rem; }}
    th, td {{ border-bottom: 1px solid #e5e5e5; padding: 0.4rem 0.5rem; text-align: left; }}
    th {{ background: #f5f7fb; }}
    .contrast-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 1rem; }}
  </style>
</head>
<body>
  <div class="hero">
    <h1>RNA-seq Analysis Portal</h1>
    <p>Output root: <code>{output_dir}</code></p>
    <p>Use this page as the delivery entrypoint for QC, differential expression, enrichment, and diagnostics.</p>
  </div>

  <div class="section">
    <h2>Quick Summary</h2>
    {preview_table(contrast_summary_path, sep=",", n=10)}
  </div>

  <div class="section">
    <h2>Upstream Provenance</h2>
    {preview_table(upstream_path, sep="\\t", n=20)}
  </div>

  <div class="section">
    <h2>Core Navigation</h2>
  <div class="grid">
    <div class="card">
      <h2>Validation</h2>
      <ul>
        <li><a href="../00_Validation/validated_inputs.txt">Validated Inputs</a></li>
        <li><a href="../00_Validation/group_summary.tsv">Group Summary</a></li>
        <li><a href="../00_Validation/contrast_summary.tsv">Contrast Summary</a></li>
        <li><a href="../00_Validation/upstream_provenance.tsv">Upstream Provenance</a></li>
      </ul>
    </div>

    <div class="card">
      <h2>QC</h2>
      <ul>
        <li><a href="../01_QC/PCA_Plot.png">PCA Plot</a></li>
        <li><a href="../01_QC/Adjusted_PCA_Plot.png">Adjusted PCA Plot</a></li>
        <li><a href="../01_QC/Sample_Distance.png">Sample Distance Heatmap</a></li>
        <li><a href="../01_QC/QC_Metadata_Associations.tsv">QC Metadata Associations</a></li>
        <li><a href="QC_Adjustment_Comparison.tsv">QC Adjustment Comparison</a></li>
      </ul>
      <a href="../01_QC/PCA_Plot.png"><img class="thumb" src="../01_QC/PCA_Plot.png" alt="PCA Plot"></a>
    </div>

    <div class="card">
      <h2>Differential Expression</h2>
      <ul>
        <li><a href="../02_DESeq2_Stats/_contrast_summary.csv">Contrast Summary</a></li>
        <li><a href="../02_DESeq2_Stats/_cooks_gene_report.tsv">Cook's Gene Report</a></li>
        <li><a href="../02_DESeq2_Stats/_cooks_sample_report.tsv">Cook's Sample Report</a></li>
        <li><a href="DEG_Counts_Barplot.png">DEG Counts Barplot</a></li>
        <li><a href="Top_DE_Genes_Heatmap.png">Top DE Genes Heatmap</a></li>
      </ul>
      <a href="DEG_Counts_Barplot.png"><img class="thumb" src="DEG_Counts_Barplot.png" alt="DEG Counts"></a>
    </div>

    <div class="card">
      <h2>Summaries</h2>
      <ul>
        <li><a href="Analysis_Summary.md">Analysis Summary</a></li>
        <li><a href="Annotated_Top_Variable_Clustermap.png">Annotated Clustermap</a></li>
        <li><a href="Sample_Outlier_Report.tsv">Sample Outlier Report</a></li>
        <li><a href="GSEA_Summary.tsv">GSEA Summary</a></li>
        <li><a href="GSEA_Summary_Dotplot.png">GSEA Summary Dotplot</a></li>
        <li><a href="GSEA_Term_Clustermap.png">GSEA Term Clustermap</a></li>
      </ul>
      <a href="GSEA_Summary_Dotplot.png"><img class="thumb" src="GSEA_Summary_Dotplot.png" alt="GSEA Dotplot"></a>
    </div>
  </div>

  <div class="section">
    <h2>Contrast Gallery</h2>
    <div class="contrast-grid">
      {''.join(contrast_cards) if contrast_cards else '<p>No contrast images found.</p>'}
    </div>
  </div>

  <div class="section">
    <h2>Outlier Preview</h2>
    {preview_table(outlier_path, sep="\\t", n=10)}
  </div>

  <div class="section">
    <h2>GSEA Preview</h2>
    {preview_table(gsea_summary_path, sep="\\t", n=12)}
  </div>
</body>
</html>
"""
    with open(os.path.join(out_dir, "Report_Index.html"), 'w', encoding='utf-8') as handle:
        handle.write(html)


def create_contrast_report_pages(output_dir, out_dir, results_dict):
    os.makedirs(out_dir, exist_ok=True)

    def rel(path):
        return os.path.relpath(path, out_dir).replace(os.sep, "/")

    for name, df in results_dict.items():
        volcano_path = os.path.join(output_dir, "03_Volcano_Plots", f"Volcano_{name}.png")
        ma_path = os.path.join(output_dir, "03_Volcano_Plots", f"MA_{name}.png")
        stat_path = os.path.join(output_dir, "02_DESeq2_Stats", f"{name}.csv")
        gsea_dirs = sorted(glob.glob(os.path.join(output_dir, "04_GSEA", name, "*")))
        motif_dir = os.path.join(output_dir, "06_Motif", name)

        clean = df.dropna(subset=['padj', 'log2FoldChange']).copy()
        top = clean.assign(
            rank_score=(-np.log10(clean['padj'].clip(lower=1e-300))) * clean['log2FoldChange'].abs()
        ).sort_values(['rank_score', 'padj'], ascending=[False, True]).head(20)

        table_html = "<p>No genes available.</p>"
        if not top.empty:
            cols = [col for col in ['log2FoldChange', 'padj', 'pvalue', 'stat'] if col in top.columns]
            headers = "".join(f"<th>{col}</th>" for col in ['gene'] + cols)
            rows = []
            for gene, row in top[cols].iterrows():
                cells = "".join(f"<td>{row[col]}</td>" for col in cols)
                rows.append(f"<tr><td>{gene}</td>{cells}</tr>")
            table_html = (
                '<div class="table-wrap"><table><thead><tr>'
                + headers
                + "</tr></thead><tbody>"
                + "".join(rows)
                + "</tbody></table></div>"
            )

        gsea_links = "".join(
            f'<li><a href="{rel(path)}">{os.path.basename(path)}</a></li>' for path in gsea_dirs
        ) or "<li>No GSEA directory found.</li>"
        motif_html = (
            f'<p><a href="{rel(motif_dir)}">Open motif results</a></p>'
            if os.path.exists(motif_dir)
            else "<p>No motif results found.</p>"
        )

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{name} - RNA-seq Contrast Report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; line-height: 1.5; color: #222; background: #fafafa; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 1rem; }}
    .card {{ border: 1px solid #ddd; border-radius: 10px; padding: 1rem; background: #fff; }}
    .thumb {{ width: 100%; border: 1px solid #eee; border-radius: 6px; margin-top: 0.75rem; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.92rem; }}
    th, td {{ border-bottom: 1px solid #e5e5e5; padding: 0.4rem 0.5rem; text-align: left; }}
    th {{ background: #f5f7fb; }}
    .table-wrap {{ overflow-x: auto; }}
    a {{ color: #0b57d0; text-decoration: none; }}
  </style>
</head>
<body>
  <p><a href="Report_Index.html">Back to portal</a></p>
  <h1>{name}</h1>
  <p><a href="{rel(stat_path)}">DE result table</a></p>

  <div class="grid">
    <div class="card">
      <h2>Volcano Plot</h2>
      <a href="{rel(volcano_path)}"><img class="thumb" src="{rel(volcano_path)}" alt="Volcano {name}"></a>
    </div>
    <div class="card">
      <h2>MA Plot</h2>
      <a href="{rel(ma_path)}"><img class="thumb" src="{rel(ma_path)}" alt="MA {name}"></a>
    </div>
  </div>

  <div class="card" style="margin-top: 1rem;">
    <h2>Top Ranked Genes</h2>
    {table_html}
  </div>

  <div class="grid" style="margin-top: 1rem;">
    <div class="card">
      <h2>GSEA Outputs</h2>
      <ul>{gsea_links}</ul>
    </div>
    <div class="card">
      <h2>Motif Outputs</h2>
      {motif_html}
    </div>
  </div>
</body>
</html>
"""
        with open(os.path.join(out_dir, f"contrast_{name}.html"), 'w', encoding='utf-8') as handle:
            handle.write(html)

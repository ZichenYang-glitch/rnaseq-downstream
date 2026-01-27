#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ================= USER CONFIGURATION =================
# Default input files (assuming you run this from the project dir)
COUNTS_FILE = 'salmon.merged.gene_counts.tsv'
METADATA_FILE = 'metadata.txt'
OUT_DIR = 'deseq2_integrated_results'

# Experimental Design Configuration
DESIGN_FACTOR = 'group'  # The column name in metadata.txt defining your groups
REFERENCE_LEVEL = 'Naive' # The baseline group (control)

# Custom Contrasts List
# Format: [('Treatment', 'Control'), ...]
# If empty [], it will default to comparing all groups against REFERENCE_LEVEL
CUSTOM_CONTRASTS = [
    # Lineage comparisons vs Naive
    ('Th1_24h', 'Naive'),
    ('Th1_72h', 'Naive'),
    ('Th2_24h', 'Naive'),
    ('Th2_72h', 'Naive'),
    # Direct lineage vs lineage
    ('Th1_24h', 'Th2_24h'),
    ('Th1_72h', 'Th2_72h'),
    # Time effect
    ('Th1_72h', 'Th1_24h'),
    ('Th2_72h', 'Th2_24h')
]

# Plotting Settings
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_context("talk")
os.makedirs(OUT_DIR, exist_ok=True)

# ================= 1. DATA LOADING =================
def load_data():
    print(f"[1/5] Loading Metadata ({METADATA_FILE}) and Counts ({COUNTS_FILE})...")
    
    # Metadata
    meta = pd.read_csv(METADATA_FILE, sep='\t')
    if 'sample_id' in meta.columns:
        meta = meta.set_index('sample_id')
    else:
        meta = meta.set_index(meta.columns[0])
    
    # Counts
    counts = pd.read_csv(COUNTS_FILE, sep='\t')
    idx_col = 'gene_name' if 'gene_name' in counts.columns else 'gene_id'
    counts = counts.groupby(idx_col).sum(numeric_only=True)
    
    # Alignment
    common_samples = [s for s in meta.index if s in counts.columns]
    meta = meta.loc[common_samples]
    counts = counts[common_samples].round().astype(int)
    
    # Filter low expressed genes
    counts = counts[counts.sum(axis=1) >= 10]
    
    return counts.T, meta

counts_T, metadata = load_data()
print(f"  Processed {counts_T.shape[0]} samples and {counts_T.shape[1]} genes.")

# ================= 2. QUALITY CONTROL =================
def run_qc(counts_T, metadata):
    print(f"[2/5] Running QC (PCA & Correlation)...")
    log_counts = np.log1p(counts_T)
    
    # A. PCA
    pca = PCA(n_components=2)
    scaled_data = StandardScaler().fit_transform(log_counts)
    pca_res = pca.fit_transform(scaled_data)
    pca_df = pd.DataFrame(pca_res, index=metadata.index, columns=['PC1', 'PC2'])
    pca_df = pd.concat([pca_df, metadata], axis=1)
    
    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=DESIGN_FACTOR, s=200, palette='Set1')
    plt.title(f"PCA: Exp Var {pca.explained_variance_ratio_[0]:.1%} / {pca.explained_variance_ratio_[1]:.1%}")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/QC_PCA_Plot.png", dpi=300)
    plt.close()
    
    # B. Correlation
    corr = log_counts.T.corr()
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr, annot=True, cmap='Blues', fmt=".2f")
    plt.title("Sample Correlation Matrix")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/QC_Correlation_Heatmap.png", dpi=300)
    plt.close()

run_qc(counts_T, metadata)

# ================= 3. DESEQ2 PROCESSING =================
print(f"[3/5] Fitting PyDESeq2 Model (Reference: {REFERENCE_LEVEL})...")
dds = DeseqDataSet(
    counts=counts_T,
    metadata=metadata,
    design_factors=DESIGN_FACTOR,
    ref_level=[DESIGN_FACTOR, REFERENCE_LEVEL],
    n_cpus=4,
    quiet=True
)
dds.deseq2()

# ================= 4. CONTRASTS & VISUALS =================
def plot_volcano(res, title, out_path):
    # Filter out genes with NA in padj or log2FoldChange
    res = res.dropna(subset=['padj', 'log2FoldChange'])
    res = res.copy() # Avoid warning
    res['sig'] = 'NS'
    res.loc[(res['padj'] < 0.05) & (res['log2FoldChange'] > 1), 'sig'] = 'Up'
    res.loc[(res['padj'] < 0.05) & (res['log2FoldChange'] < -1), 'sig'] = 'Down'
    
    plt.figure(figsize=(8, 6))
    colors = {'NS': 'lightgrey', 'Up': '#e41a1c', 'Down': '#377eb8'}
    sns.scatterplot(data=res, x='log2FoldChange', y=-np.log10(res['padj']), 
                    hue='sig', palette=colors, alpha=0.6, s=15, edgecolor=None)
    
    plt.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    plt.axvline(1, color='black', linestyle='--', alpha=0.5)
    plt.axvline(-1, color='black', linestyle='--', alpha=0.5)
    
    n_up = sum(res['sig'] == 'Up')
    n_down = sum(res['sig'] == 'Down')
    plt.title(f"{title}\n(Up: {n_up}, Down: {n_down})")
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(padj)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

print(f"[4/5] Running Contrasts...")
summary_list = []

# If no custom contrasts, compare everything vs reference
if not CUSTOM_CONTRASTS:
    all_groups = metadata[DESIGN_FACTOR].unique()
    target_contrasts = [(g, REFERENCE_LEVEL) for g in all_groups if g != REFERENCE_LEVEL]
else:
    target_contrasts = CUSTOM_CONTRASTS

for treat, ctrl in target_contrasts:
    print(f"  -> {treat} vs {ctrl}")
    try:
        stat_res = DeseqStats(dds, contrast=[DESIGN_FACTOR, treat, ctrl], n_cpus=4, quiet=True)
        stat_res.summary()
        res_df = stat_res.results_df
        
        # Save CSV
        name = f"{treat}_vs_{ctrl}"
        res_df.to_csv(f"{OUT_DIR}/{name}.csv")
        
        # Plot Volcano
        plot_volcano(res_df, name, f"{OUT_DIR}/Volcano_{name}.png")
        
        # Stats for summary
        summary_list.append({
            'Contrast': name,
            'Up': sum((res_df['padj'] < 0.05) & (res_df['log2FoldChange'] > 1)),
            'Down': sum((res_df['padj'] < 0.05) & (res_df['log2FoldChange'] < -1))
        })
    except Exception as e:
        print(f"     Error running contrast {treat} vs {ctrl}: {e}")

# ================= 5. FINAL REPORT =================
print(f"[5/5] Finishing Up...")
summary_df = pd.DataFrame(summary_list)
summary_df.to_csv(f"{OUT_DIR}/_summary_report.csv", index=False)
print("\nSummary of results:")
print(summary_df)
print(f"\nDone! Results are in: {os.path.abspath(OUT_DIR)}")

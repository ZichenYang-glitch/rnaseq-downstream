import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def run_qc(counts_T, metadata, design_col, out_dir):
    """Generates PCA and Correlation plots."""
    os.makedirs(out_dir, exist_ok=True)
    log_counts = np.log1p(counts_T)
    
    # PCA
    pca = PCA(n_components=2)
    scaled = StandardScaler().fit_transform(log_counts)
    coords = pca.fit_transform(scaled)
    pca_df = pd.DataFrame(coords, index=metadata.index, columns=['PC1', 'PC2'])
    pca_df = pd.concat([pca_df, metadata], axis=1)
    
    plt.figure(figsize=(8,6))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=design_col, style=design_col, s=150)
    plt.title(f"PCA Plot")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'PCA_Plot.png'), dpi=300)
    plt.close()
    
    # Correlation
    corr = log_counts.T.corr()
    plt.figure(figsize=(10,8))
    sns.heatmap(corr, cmap='RdBu_r', center=0.8) # Red-Blue color map
    plt.title("Sample Correlation")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'Sample_Correlation.png'), dpi=300)
    plt.close()

def fit_model(counts_T, metadata, design_col, ref_level):
    print("  Fitting DESeq2 model...")
    dds = DeseqDataSet(
        counts=counts_T, 
        metadata=metadata, 
        design_factors=design_col,
        ref_level=[design_col, ref_level],
        n_cpus=4, quiet=True
    )
    dds.deseq2()
    return dds

def run_contrasts(dds, contrasts, design_col, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    results = {}
    
    for treat, ctrl in contrasts:
        print(f"  Comparing: {treat} vs {ctrl}")
        try:
            stat_res = DeseqStats(dds, contrast=[design_col, treat, ctrl], n_cpus=4, quiet=True)
            stat_res.summary()
            df = stat_res.results_df
            df.to_csv(os.path.join(out_dir, f"{treat}_vs_{ctrl}.csv"))
            results[f"{treat}_vs_{ctrl}"] = df
        except Exception as e:
            print(f"    [Warning] Failed contrast {treat} vs {ctrl}: {e}")
            
    return results

def plot_volcano(results, out_dir, padj_thresh=0.05, logfc_thresh=1.0):
    os.makedirs(out_dir, exist_ok=True)
    for name, df in results.items():
        df = df.dropna(subset=['padj', 'log2FoldChange'])
        cols = ['grey'] * len(df)
        
        # Vectorized color assignment
        # Red: Sig Up
        mask_up = (df['padj'] < padj_thresh) & (df['log2FoldChange'] > logfc_thresh)
        # Blue: Sig Down
        mask_down = (df['padj'] < padj_thresh) & (df['log2FoldChange'] < -logfc_thresh)
        
        plt.figure(figsize=(6,5))
        plt.scatter(df['log2FoldChange'], -np.log10(df['padj']), c='lightgrey', s=5, alpha=0.5)
        plt.scatter(df.loc[mask_up, 'log2FoldChange'], -np.log10(df.loc[mask_up, 'padj']), c='red', s=10, label='Up')
        plt.scatter(df.loc[mask_down, 'log2FoldChange'], -np.log10(df.loc[mask_down, 'padj']), c='blue', s=10, label='Down')
        
        plt.axhline(-np.log10(padj_thresh), ls='--', c='k', lw=0.5)
        plt.axvline(logfc_thresh, ls='--', c='k', lw=0.5)
        plt.axvline(-logfc_thresh, ls='--', c='k', lw=0.5)
        
        plt.title(name)
        plt.xlabel("log2 Fold Change")
        plt.ylabel("-log10(padj)")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"Volcano_{name}.png"), dpi=300)
        plt.close()

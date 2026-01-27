import pandas as pd
import os

def create_master_table(metadata, design_col, results_dict, tpm_file, counts_T, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Base Expression (TPM or Counts)
    if os.path.exists(tpm_file):
        base = pd.read_csv(tpm_file, sep='\t')
        idx = 'gene_name' if 'gene_name' in base.columns else 'gene_id'
        if idx not in base.columns: idx = base.columns[0]
        base = base.groupby(idx).mean(numeric_only=True)
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
    for name, df in results_dict.items():
        sub = df[['log2FoldChange', 'padj']].copy()
        sub.columns = [f"log2FC_{name}", f"padj_{name}"]
        final = final.join(sub, how='left')
        
    out_path = os.path.join(out_dir, "Master_Expression_Table.csv")
    final.to_csv(out_path)
    print(f"  Master table saved to: {out_path}")


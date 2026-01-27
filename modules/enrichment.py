import gseapy as gp
import pandas as pd
import os

def run_gsea(results_dict, gene_sets, out_dir):
    """Runs GSEA Preranked for each contrast."""
    os.makedirs(out_dir, exist_ok=True)
    
    for name, df in results_dict.items():
        print(f"  Running GSEA for: {name}")
        
        # Prepare rank
        rk = df[['log2FoldChange']].copy()
        rk['gene'] = df.index
        rk = rk.dropna().sort_values('log2FoldChange', ascending=False)
        rk = rk[['gene', 'log2FoldChange']]
        
        if len(rk) < 100:
            print("    [Skip] Not enough genes.")
            continue
            
        for gs in gene_sets:
            try:
                gp.prerank(rnk=rk, 
                           gene_sets=gs, 
                           outdir=os.path.join(out_dir, name, gs),
                           format='png', 
                           seed=42, 
                           verbose=False)
            except Exception as e:
                print(f"    [Error] {gs}: {e}")

import os

import gseapy as gp


def run_gsea(results_dict, gene_sets, out_dir, permutations=1000, rank_metric='stat'):
    """Runs GSEA Preranked for each contrast."""
    os.makedirs(out_dir, exist_ok=True)
    
    for name, df in results_dict.items():
        print(f"  Running GSEA for: {name}")
        
        # Prepare rank
        metric = rank_metric if rank_metric in df.columns else 'log2FoldChange'
        rk = df[[metric]].copy()
        rk['gene'] = df.index
        rk = rk.dropna().sort_values(metric, ascending=False)
        rk = rk[['gene', metric]]
        
        if len(rk) < 100:
            print("    [Skip] Not enough genes.")
            continue
            
        for gs in gene_sets:
            try:
                gp.prerank(rnk=rk, 
                           gene_sets=gs, 
                           outdir=os.path.join(out_dir, name, gs),
                           permutation_num=permutations,
                           format='png', 
                           seed=42, 
                           verbose=False)
            except Exception as e:
                print(f"    [Error] {gs}: {e}")

import os
import subprocess
import shutil
import pandas as pd

def check_homer_installed():
    """Checks if findMotifs.pl is in the system PATH."""
    return shutil.which('findMotifs.pl') is not None

def run_homer_motif(gene_list_file, out_dir, species='human', region_args=['-start', '-400', '-end', '100']):
    """
    Runs HOMER findMotifs.pl on a gene list.
    
    Args:
        gene_list_file (str): Path to text file with gene IDs (Entrez or Symbol).
        out_dir (str): Directory to save HOMER results.
        species (str): 'human', 'mouse', etc.
        region_args (list): List of arguments for promoter region (e.g. ['-start', '-400']).
    """
    if not check_homer_installed():
        print("  [Warning] HOMER (findMotifs.pl) not found. Skipping.")
        return

    os.makedirs(out_dir, exist_ok=True)
    
    # Construct command
    # findMotifs.pl <inputList.txt> <species> <outputDir> [options]
    cmd = ['findMotifs.pl', gene_list_file, species, out_dir] + region_args + ['-p', '4']
    
    print(f"  Running HOMER: {' '.join(cmd)}")
    
    try:
        # Run command, redirecting stdout/stderr to a log file inside the output dir
        with open(os.path.join(out_dir, 'homer.log'), 'w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as e:
        print(f"  [Error] HOMER failed for {gene_list_file}. See log in {out_dir}.")

def run_motif_analysis(results, out_dir, species, padj_thresh=0.05, logfc_thresh=1.0):
    """
    Main function to orchestrate Motif analysis for all contrasts.
    Splits DEGs into UP and DOWN lists and runs HOMER for each.
    """
    if not check_homer_installed():
        print("[Warning] HOMER is not installed or not in PATH. Skipping Motif analysis.")
        return

    print(f"[Motif] Starting HOMER analysis for species: {species}")
    
    for name, df in results.items():
        # Prepare sub-directories
        contrast_dir = os.path.join(out_dir, name)
        
        # Filter significant genes
        # We need clean gene names (index of df). 
        # HOMER can handle Gene Symbols or Entrez IDs.
        
        # UP genes
        up_genes = df[(df['padj'] < padj_thresh) & (df['log2FoldChange'] > logfc_thresh)].index.tolist()
        if len(up_genes) > 5: # Only run if enough genes
            up_dir = os.path.join(contrast_dir, 'UP')
            os.makedirs(up_dir, exist_ok=True)
            
            # Save gene list
            up_list_file = os.path.join(up_dir, 'genes_up.txt')
            with open(up_list_file, 'w') as f:
                f.write('\n'.join(up_genes))
                
            print(f"  Processing {name} UP ({len(up_genes)} genes)...")
            run_homer_motif(up_list_file, up_dir, species)
        
        # DOWN genes
        down_genes = df[(df['padj'] < padj_thresh) & (df['log2FoldChange'] < -logfc_thresh)].index.tolist()
        if len(down_genes) > 5:
            down_dir = os.path.join(contrast_dir, 'DOWN')
            os.makedirs(down_dir, exist_ok=True)
            
            # Save gene list
            down_list_file = os.path.join(down_dir, 'genes_down.txt')
            with open(down_list_file, 'w') as f:
                f.write('\n'.join(down_genes))
                
            print(f"  Processing {name} DOWN ({len(down_genes)} genes)...")
            run_homer_motif(down_list_file, down_dir, species)

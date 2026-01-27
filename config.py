
# ================= CONFIGURATION =================

# 1. Input File Paths (Relative to project root or absolute)
COUNTS_FILE = 'salmon.merged.gene_counts.tsv'
TPM_FILE    = 'salmon.merged.gene_tpm.tsv'
METADATA_FILE = 'metadata.txt'

# 2. Output Settings
OUTPUT_DIR = 'Final_Analysis_Pipeline'

# 3. Experimental Design
# Column in metadata.txt
DESIGN_FACTOR = 'group' 
# Control group
REFERENCE_LEVEL = 'M0'

# 4. Contrasts [('Treatment', 'Control'), ...]
CONTRASTS = [
    # Vs M0
    ('M1', 'M0'),
    ('M2', 'M0'),
    
    
    # Direct Comparisons
    ('M1', 'M2'),
    
]

# 5. Analysis Parameters
MIN_COUNTS = 10        # Filter genes with sum < 10
PADJ_THRESH = 0.05     # Significance threshold
LOGFC_THRESH = 1.0     # Log2FC threshold for coloring

# 6. GSEA Settings
RUN_GSEA = True
GSEA_GENE_SETS = ['KEGG_2019_Mouse', 'GO_Biological_Process_2021']
GSEA_PERMUTATIONS = 1000

# 7. Motif Analysis Settings (HOMER)
RUN_MOTIF = True
HOMER_SPECIES = 'mouse' # 'human', 'mouse', etc.

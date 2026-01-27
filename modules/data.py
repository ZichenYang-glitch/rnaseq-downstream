import pandas as pd
import sys

def load_metadata(path, design_col):
    """Loads and validates metadata."""
    try:
        df = pd.read_csv(path, sep='\t')
        # Handle index
        if 'sample_id' in df.columns:
            df = df.set_index('sample_id')
        else:
            df = df.set_index(df.columns[0])
            
        if design_col not in df.columns:
            raise ValueError(f"Design column '{design_col}' not found in metadata.")
            
        return df
    except Exception as e:
        print(f"[Error] Loading metadata: {e}")
        sys.exit(1)

def load_counts(path, metadata_samples, min_counts=10):
    """Loads counts, aligns with metadata, and filters."""
    try:
        df = pd.read_csv(path, sep='\t')
        # Identify gene column
        idx = 'gene_name' if 'gene_name' in df.columns else 'gene_id'
        if idx not in df.columns: idx = df.columns[0]
        
        df = df.groupby(idx).sum(numeric_only=True)
        
        # Align samples
        common = [s for s in metadata_samples if s in df.columns]
        if not common:
            raise ValueError("No matching samples between metadata and counts!")
            
        df = df[common].fillna(0).round().astype(int)
        
        # Filter
        df = df[df.sum(axis=1) >= min_counts]
        
        return df.T  # Return samples x genes
    except Exception as e:
        print(f"[Error] Loading counts: {e}")
        sys.exit(1)

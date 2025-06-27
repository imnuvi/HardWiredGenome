import pandas as pd
import numpy as np
import anndata
from mapping import generate_gene_id_name_map

def fetch_genexpression_data(a_matrix_adata, celltype, cell):
    gene_id_name_map, gene_name_id_map = generate_gene_id_name_map()
    gex_adata = anndata.read_h5ad('./data/Genexpression.h5ad')
    marker_gene_df = pd.read_csv('./data/clean_marker_genes.csv')
    
    fibroblast_markers = marker_gene_df[marker_gene_df['cell_type'] == celltype]
    marker_genes = list(fibroblast_markers['gene_name'])
    
    
    marker_ids = [ gene_name_id_map[i] for i in marker_genes ]
    expr_values = gex_adata.X[gex_adata.obs_names == cell].flatten()
    
    expr_series = pd.Series(expr_values, index=gex_adata.var_names)
    marker_expr = expr_series[expr_series.index.isin(marker_ids)]
    
    # threshold = marker_expr.mean()
    threshold = marker_expr.median()
    high_expr = expr_series[expr_series > threshold]
    top_genes = high_expr.sort_values(ascending=False)
    top_genes = top_genes[top_genes.index.notna()]
    fibroblast_expressed_genes = top_genes.index.tolist()
    
    fibroblast_adata = a_matrix_adata[a_matrix_adata.obs_names.isin(fibroblast_expressed_genes), a_matrix_adata.var_names.isin(fibroblast_expressed_genes)].copy()
    return fibroblast_adata

fetch_genexpression_data(a_matrix_adata, 'fibroblasts', 'Fibroblast'
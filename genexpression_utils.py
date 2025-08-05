import pandas as pd
import numpy as np
import anndata
import scanpy as sc

from mapping import generate_gene_id_name_map, get_a_matrix_threshold

def fetch_genexpression_data(a_matrix_adata, gex_adata, marker_genes, cell):
    gene_id_name_map, gene_name_id_map = generate_gene_id_name_map()

    marker_ids = [ ]
    for i in marker_genes:
        try:
            gid = gene_name_id_map[i]
            marker_ids.append(gid)
        except:
            print(f'no mapping found for {i}')
            continue
        
    expr_values = gex_adata.X[gex_adata.obs_names == cell].flatten()
    
    expr_series = pd.Series(expr_values, index=gex_adata.var_names)
    marker_expr = expr_series[expr_series.index.isin(marker_ids)]

    threshold = marker_expr.mean()
    # threshold = marker_expr.median()
    high_expr = expr_series[expr_series > threshold]
    top_genes = high_expr.sort_values(ascending=False)
    top_genes = top_genes[top_genes.index.notna()]
    fibroblast_expressed_genes = top_genes.index.tolist()
    
    fibroblast_adata = a_matrix_adata[a_matrix_adata.obs_names.isin(fibroblast_expressed_genes), a_matrix_adata.var_names.isin(fibroblast_expressed_genes)].copy()
    return fibroblast_adata

def fetch_fibroblast_data():
    a_matrix_adata, b_matrix_true = get_a_matrix_threshold(300)
    
    gex_adata = anndata.read_h5ad('./data/Genexpression.h5ad')
    marker_gene_df = pd.read_csv('./data/clean_marker_genes.csv')

    celltype = 'fibroblasts'
    obs_name = 'Fibroblast'
    
    fibroblast_markers = marker_gene_df[marker_gene_df['cell_type'] == celltype]
    marker_genes = list(fibroblast_markers['gene_name'])

    fib_adata = fetch_genexpression_data(a_matrix_adata, gex_adata, marker_genes, obs_name)
    return fib_adata

def fetch_myotube_data():
    a_matrix_adata, b_matrix_true = get_a_matrix_threshold(300)
    gex_adata = anndata.read_h5ad('./data/Genexpression.h5ad')
    
    marker_gene_sets = get_marker_gene_sets()
    myotube_markers = set(marker_gene_sets['Muscle'].tolist())
    
    marker_genes = list(myotube_markers)
    if 'nan' in marker_genes:
        marker_genes.remove('nan')
    
    print(myotube_markers)
    obs_name = 'Myotube'
    muscle_adata = fetch_genexpression_data(a_matrix_adata, gex_adata, marker_genes, obs_name)
    return muscle_adata

def fetch_ipsc_data():
    a_matrix_adata, b_matrix_true = get_a_matrix_threshold(300)
    gex_adata = anndata.read_h5ad('./data/Genexpression.h5ad')
    
    marker_gene_sets = get_marker_gene_sets()
    myotube_markers = set(marker_gene_sets['Muscle'].tolist())
    
    marker_genes = list(myotube_markers)
    if 'nan' in marker_genes:
        marker_genes.remove('nan')
    
    print(myotube_markers)
    obs_name = 'Myotube'
    muscle_adata = fetch_genexpression_data(a_matrix_adata, gex_adata, marker_genes, obs_name)
    return muscle_adata

def get_marker_gene_sets():
    return pd.read_csv('./data/gene_sets.tsv', sep='\t')

def process_gex_anndata(fibroblast, fib_marker_ids, thresh):
    a_matrix_adata, b_matrix_true = get_a_matrix_threshold(300)

    sc.pp.normalize_total(fibroblast, target_sum=1e4)
    sc.pp.log1p(fibroblast)
    
    fibroblast.var['GeneName'] = fibroblast.var_names
    fibroblast.var_names = fibroblast.var['gene_id']

    gex_adata = fibroblast

    
    
    expr_values = gex_adata.X.flatten()
    
    expr_series = pd.Series(expr_values, index=gex_adata.var_names)

    if fib_marker_ids:
        marker_expr = expr_series[expr_series.index.isin(fib_marker_ids)]

    if thresh:
        threshold = thresh
    else:
        threshold = marker_expr.mean()
    
    print(threshold)
    high_expr = expr_series[expr_series > threshold]
    top_genes = high_expr.sort_values(ascending=False)
    top_genes = top_genes[top_genes.index.notna()]
    fibroblast_expressed_genes = top_genes.index.tolist()
    
    fibroblast_adata = a_matrix_adata[a_matrix_adata.obs_names.isin(fibroblast_expressed_genes), a_matrix_adata.var_names.isin(fibroblast_expressed_genes)].copy()
    
    return fibroblast_adata
    
def fetch_hyb_nodes():
    gene_id_name_map, gene_name_id_map = generate_gene_id_name_map()

    
    HYB_anndata = fetch_HYB_genexpression()

    fibroblast = HYB_anndata[HYB_anndata.obs_names == 'Control']
    reprogrammed = HYB_anndata[HYB_anndata.obs_names == 'mmMYOD1']
    

    marker_gene_sets = get_marker_gene_sets()
    myotube_markers = set(marker_gene_sets['Muscle'].tolist())
    
    muscle_marker_genes = list(myotube_markers)
    
    
    fib_markers = set(marker_gene_sets['Fibroblast'].tolist())
    
    fib_marker_genes = list(fib_markers)
    
    fib_marker_ids = [ ]
    for i in fib_marker_genes:
        try:
            gid = gene_name_id_map[i]
            fib_marker_ids.append(gid)
        except:
            print(f'no mapping found for {i}')
            continue
    
    muscle_marker_ids = [ ]
    for i in fib_marker_genes:
        try:
            gid = gene_name_id_map[i]
            muscle_marker_ids.append(gid)
        except:
            print(f'no mapping found for {i}')
            continue
    thresh = 0.38259256
    hyb_fib = process_gex_anndata(fibroblast, fib_marker_ids, None)
    hyb_reprogrammed = process_gex_anndata(reprogrammed, None, thresh)
    myotube = fetch_myotube_data()

    return hyb_fib, hyb_reprogrammed, myotube
    

def fetch_HYB_genexpression():
    hybrid_expression_path = '/scratch/indikar_root/indikar1/shared_data/HYB/anndata/pseudobulk.h5ad'
    hybrid_adata = anndata.read_h5ad(hybrid_expression_path)
    return hybrid_adata

if __name__ == '__main__':
    # fib_adata = fetch_fibroblast_data()
    # print(fib_adata)
    myotube_adata = fetch_myotube_data()
    print(myotube_adata)
    
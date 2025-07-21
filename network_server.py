import pandas as pd
import numpy as np
import anndata
from scipy import sparse 
from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import json

from flask import Flask, request, render_template

from mapping import get_sorted_gene_order, generate_gene_id_name_map, get_master_regulator_list, get_TF_lists, get_a_matrix_threshold
from genexpression_utils import fetch_hyb_nodes
from constants import OPERATIONS_DIRECTORY
from network_utils import  GeneUtils


app = Flask(__name__)


GUtils = GeneUtils()
hyb_fib, hyb_reprogrammed, myotube = fetch_hyb_nodes()

hyb_fib_genes = hyb_fib.var_names
hyb_reprogrammed_genes = hyb_reprogrammed.var_names
myotube_genes = myotube.var_names

# # Load Content to memory
# a_matrix_adata, b_matrix_true = get_a_matrix_threshold(300)
# b_matrix = anndata.read_h5ad('/scratch/indikar_root/indikar1/shared_data/HWG/operations/B_matrices_true.h5ad')
# master_regulator_list = get_master_regulator_list()
# repressorlist, activatorlist, conflictedlist, tf_list = get_TF_lists()
# print(f'total Transcription Factors : {len(set(tf_list))}')
# genes = get_sorted_gene_order()
# b_matrix_true = a_matrix_adata[a_matrix_adata.obs_names.isin(tf_list)]



# # Genexpression Vectors
# gex_adata = anndata.read_h5ad('./data/Genexpression.h5ad')
# marker_gene_df = pd.read_csv('./data/clean_marker_genes.csv')

# fibroblast_markers = marker_gene_df[marker_gene_df['cell_type'] == 'fibroblasts']
# print(fibroblast_markers)
# marker_genes = list(fibroblast_markers['gene_name'])
# gene_id_name_map, gene_name_id_map = generate_gene_id_name_map()

# marker_ids = [ gene_name_id_map[i] for i in marker_genes ]
# print(len(marker_ids), len(set(marker_ids)))
# print(gex_adata.obs_names)
# expr_values = gex_adata.X[gex_adata.obs_names == 'Fibroblast'].flatten()
# print(expr_values)
# expr_series = pd.Series(expr_values, index=gex_adata.var_names)
# marker_expr = expr_series[expr_series.index.isin(marker_ids)]
# print(len(marker_expr))




# # threshold = marker_expr.mean()
# threshold = marker_expr.median()
# high_expr = expr_series[expr_series > threshold]
# top_genes = high_expr.sort_values(ascending=False)
# top_genes = top_genes[top_genes.index.notna()]
# fibroblast_expressed_genes = top_genes.index.tolist()



# fibroblast_adata = a_matrix_adata[a_matrix_adata.obs_names.isin(fibroblast_expressed_genes), a_matrix_adata.var_names.isin(fibroblast_expressed_genes)].copy()

def adata_to_json(adata):
    coo = adata.X.tocoo()
    entries = np.vstack((coo.row, coo.col)).T.tolist()
    gene_id_name_map = GUtils.gene_id_name_map
    return {
        "rowLabels": [gene_id_name_map.get(i, i) for i in adata.obs_names.tolist()],
        "colLabels": [gene_id_name_map.get(i, i) for i in adata.var_names.tolist()],
        "entries": entries,        # sparse 1-positions
        "nRows": [adata.n_obs],
        "nCols": adata.n_vars
    }

@app.route("/fetch_matrix")
def fetch_matrix():
    adata = GUtils.a_matrix_adata
    return adata_to_json(adata)

@app.route("/matrix", methods=['GET'])
def matrix():
    return render_template('matrix_viewer.html')


@app.route("/fetch_network", methods=['GET', 'POST'])
def index():
    startindexx = int(request.args['start_index_x'])
    startindexy = int(request.args['start_index_y'])
    endindexx = int(request.args['end_index_x'])
    endindexy = int(request.args['end_index_y'])
    perturbations = list(request.args.get('perturbations','').split(';'))
    celltype = request.args.get('celltype','')



    ntype = 'base'
    genes_subx = GUtils.genes[startindexx:endindexx]
    genes_suby = GUtils.genes[startindexy:endindexy]
    hwg_subnetx = []
    hwg_subnety = []

    if celltype == 'fibroblast':

        print('within fibroblast')
        print(len(genes_subx), len(genes_suby))
        hwg_subnetx = genes_subx
        hwg_subnety = genes_suby
        genes_subx = [gene for gene in genes_subx if gene in hyb_fib_genes]
        genes_suby = [gene for gene in genes_suby if gene in hyb_fib_genes]
        print(len(genes_subx), len(genes_suby))




    # full_network_G = GUtils.construct_network(genes_sub, [])
    full_network_G = GUtils.construct_network_x_y(genes_subx, genes_suby, hwg_subnetx, hwg_subnety)

    full_network_G = GUtils.add_network_props(full_network_G, ntype)

    if len(perturbations) > 0 and perturbations[0] != '':
        full_network_G = GUtils.insert_perturbation(full_network_G, perturbations, GUtils.a_matrix_adata)

    network_res = GUtils.export_graph_as_json(full_network_G)
    network_json = json.dumps(network_res)

    return network_json
    
@app.route("/network", methods=['GET', 'POST'])
def network():
    startindexx = int(request.args['start_index_x'])
    endindexx = int(request.args['end_index_x'])
    startindexy = int(request.args['start_index_y'])
    endindexy = int(request.args['end_index_y'])
    perturbations = request.args.get('perturbations', '')
    celltype = request.args.get('celltype', '')

    # startindexx = 1000
    # startindexy = 2000
    # endindexx = 1000
    # endindexy = 2000
    print('-------------------------------')
    print(startindexx, startindexy)
    return render_template('graph_viewer.html', startindexx=startindexx, startindexy=startindexy, endindexx=endindexx, endindexy=endindexy, perturbations=perturbations, celltype=celltype)




@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"





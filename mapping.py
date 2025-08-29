from constants import DATA_DIRECTORY,  STRING_ALIAS_PATH, STRING_URL_PATH, ENSEMBL_PROTEIN_GENE_PATH, STRING_EXTRACTED_ALIAS_PATH, STRING_EXTRACTED_URL_PATH, ENSEMBL_MAPPING_PATH, UNIQUE_GENE_SET, \
    HURI_URL_PATH, HI_UNION_PATH, LIT_BM_PATH, HUMAN_TF_PATH, HARD_WIRED_GENOME_A_MATRIX_PATH, STRING_PROTEIN_GENE_PATH, STRING_UNIQUE_GENE_SET, HURI_UNIQUE_GENE_SET, HUMAN_TF_SET_PATH, \
    HARD_WIRED_GENOME_B_MATRIX_PATH, HWG_BASE_PATH, ENSEMBL_BASE_PATH, STRING_PROTEIN_SET_FROM_LINKS, STRING_UNIQUE_PROTEIN_SET, modify_HWG_path, HUMAN_TF_SET_PATH, HUMAN_TF_IDLIST_PATH, TF_MOTIF_PATH
from constants import HWG_BASE_PATH, HI_UNION_PATH, LIT_BM_PATH, ENSEMBL_PROTEIN_GENE_PATH, OPERATIONS_DIRECTORY, HWG_B_MATRICES, HWG_A_MATRICES, HTF_MOTIFS_DIR, CISBP_MOTIFS_DIR

import re
import os
import numpy as np
import pandas as pd
import anndata
import scipy.sparse as sp
from anndata import AnnData
from Bio import motifs

def fetch_ensembl_mapping():
    from biomart import BiomartServer
    """ 
    Gets The Full Ensembl mapping data from biomart
    """
    atts = ['external_gene_name','external_gene_source','ensembl_gene_id',
            'ensembl_transcript_id','ensembl_peptide_id', 'uniprot_gn_id', 'uniprot_gn_symbol']

    server = BiomartServer( "http://useast.ensembl.org/biomart" )
    #server = BiomartServer("http://www.biomart.org/biomart")
    hge = server.datasets['hsapiens_gene_ensembl']

    if not os.path.exists(ENSEMBL_BASE_PATH):
        os.makedirs(ENSEMBL_BASE_PATH)

#   Uncomment to get all avaliable biomart attributes. Not really a good way to get mapping data. 
#     from pprint import pprint
#     pprint(hge.attributes)

    s = hge.search({'attributes': atts}, header=1)

    # ensembl_mapping_path = f"{DATA_DIRECTORY}/ensembl_mapping.txt"

    with open(ENSEMBL_MAPPING_PATH, 'wb') as file:
        for l in s.iter_lines():
            file.write(l + b'\n')
        print(f"Ensembl mapping data saved to {ENSEMBL_MAPPING_PATH}")

def generate_string_mappings():
    """
    We Consider ENSEMBL ids as primary identifiers and convert every other identifier to ensembl
    """

    # Takes a long time and is not really useful. Same thing is done by fetch_ensembl_mapping()
    
    # Find Unique Protein interaction identifiers in STRING
    stringfile = STRING_EXTRACTED_URL_PATH
    lineset = set()
    with open(stringfile, 'r') as f:
        count = 0
        for line in f:
            count += 1
            # print(f'running line {count}')
            p1, p2, score = line.strip().split(' ')
            lineset = lineset | {p1, p2}

    with open(STRING_PROTEIN_SET_FROM_LINKS, "w") as f:
        for item in lineset:
            f.write(f"{item}\n")
    return lineset

def generate_gene_id_name_map():
    ENSEMBL_GENE_NAME_PATH = "./data/Gene Names.tsv"

    usecols = ['Gene stable ID', 'Gene name']
    with open(ENSEMBL_GENE_NAME_PATH, 'r') as f:
        ensembl_gene_name_path = f.readlines()
    # ensembl_df = pd.read_csv(ENSEMBL_GENE_NAME_PATH, sep='\t', dtype=str, usecols=usecols, low_memory=False)

    ensmap = {}
    for row in ensembl_gene_name_path[1:]:
        splitrow = row.strip().split('\t')
        if len(splitrow) > 7:
            gene_id = splitrow[0]
            gene_name = splitrow[-1]
            ensmap[gene_id] = gene_name
                
    print(f'Extracted {len(ensmap)} mappings from {len(ensembl_gene_name_path)} ensembl lines')
        
    # gene_id_name_map = dict(zip(ensembl_df['Gene stable ID'], ensembl_df['Gene name']))
    
    gene_id_name_map_gtf = pd.read_csv('./data/gene_id_name_map_gtf.csv')
    gene_id_name_map_dict = dict(zip(gene_id_name_map_gtf['gene_id'], gene_id_name_map_gtf['gene_name']))
    print(f'Extracted {len(gene_id_name_map_dict)} mappings from gene id gtf file')
    
    gene_id_name_map_dict.update(ensmap)
    
    name_id_map = {}
    for key, value in gene_id_name_map_dict.items():
        if value:
            name_id_map[value] = key

    print(f'Total Mappings Extracted {len(gene_id_name_map_dict)}')
    
    return gene_id_name_map_dict, name_id_map
        
def generate_protein_gene_mappings():
    """
    Generates mapping of Proteins to Genes
    """
    unique_genes = set()
    with open(ENSEMBL_MAPPING_PATH, "r") as infile, open(ENSEMBL_PROTEIN_GENE_PATH, "w") as outfile:
        for line in infile:
            if "ENSP" in line and "ENSG" in line:
                matches = re.findall(r"\b(ENSP\w+|ENSG\w+)", line)
                gene = matches[0]
                unique_genes.add(gene)
                outfile.write(' '.join(matches) + '\n')

    print(f"Total unique genes from Ensembl: {len(unique_genes)}")
    with open(UNIQUE_GENE_SET, "w") as gene_file:
        for item in unique_genes:
            gene_file.write(f"{item}\n")

def generate_protein_gene_mapping_STRING():
    """
    Generates mapping of Proteins to Genes from STRING dataset
    """
    protein_gene_map = {}
    with open(STRING_EXTRACTED_ALIAS_PATH, "r") as aliasfile, open(STRING_PROTEIN_GENE_PATH, "w") as outfile:
        unique_genes = set()
        for line in aliasfile:
            if "Ensembl_gene" in line:
                matches = re.findall(r"\b(ENSP\w+|ENSG\w+)", line)
                gene = matches[1]
                protein = matches[0]
                protein_gene_map[protein] = gene
                unique_genes.add(gene)
                outfile.write(' '.join(matches) + '\n')

    print(f"Total unique genes STRING: {len(unique_genes)}")
    with open(STRING_UNIQUE_GENE_SET, "w") as gene_file:
        for item in unique_genes:
            gene_file.write(f"{item}\n")

    return protein_gene_map

        
def generate_unique_gene_set_Huri():
    with open(HI_UNION_PATH, "r") as huri_links, open(LIT_BM_PATH, "r") as litbm_links, open(HURI_UNIQUE_GENE_SET, "a") as gene_file:
        unique_genes = set()
        for link in huri_links:
            g1, g2 = link.strip().split()
            unique_genes.add(g1)
            unique_genes.add(g2)
            
        print(f"Unique Genes from HI-Union: {len(unique_genes)}")

        unique_litbm_genes = set()
        for link in litbm_links:
            g1, g2 = link.strip().split()
            unique_litbm_genes.add(g1)
            unique_litbm_genes.add(g2)
            
        print(f"Unique Genes from LIT-BM: {len(unique_litbm_genes)}")
        
        
        unique_genes = unique_genes | unique_litbm_genes
        print(f"Total unique genes HURI ( union of HI-Union + LIT-BM ): {len(unique_genes)}")
        for item in unique_genes:
            gene_file.write(f"{item}\n")
    return unique_genes

def get_unique_genes():

    string_protein_gene_map = get_string_protein_gene_map()

    with open(f'{DATA_DIRECTORY}/STRING/unique_protein_set.txt', "r") as f:
        # print("Unique Proteins from STRING:", f.readlines()[:10])
        lines =f.readlines()[1:]

        prots = [ ]
        for line in lines:
            protid = line.strip().split('.')
            if len(protid) > 1:
                prots.append(protid[1])
        proteins = [string_protein_gene_map.get(i) for i in prots]
        print(proteins[:10])
        string_nodes = set(proteins)
    print(list(string_nodes)[:5], len(string_nodes))


    with open(HURI_UNIQUE_GENE_SET, "r") as f:
        huri_nodes = set(line.strip() for line in f)

    nodes = string_nodes | huri_nodes
    print(f"Unique genes from STRING: {len(string_nodes)}")
    print(f"Unique genes from HURI: {len(huri_nodes)}")
    print("Total unique genes:", len(nodes))
    nodes = sorted(nodes)
    return nodes

def generate_human_tf_set():
    tf_df = pd.read_csv(HUMAN_TF_PATH)
    print(tf_df.head())
    print(tf_df.columns)
    subset = tf_df[["Ensembl ID", "Is TF?"]]

    print(subset.head())

    subset["Is TF?"] = subset["Is TF?"].replace({"Yes": 1, "No": 0})

    subset.rename(columns={"Is TF?": "TF", "Ensembl ID": "Gene"}, inplace=True)

    subset.to_csv(HUMAN_TF_SET_PATH, index=False)

def get_string_protein_gene_map():
    with open(STRING_PROTEIN_GENE_PATH, "r") as protein_gene_file:
        string_protein_gene_map = {}
        for line in protein_gene_file:
            protein, gene = line.split()
            string_protein_gene_map[protein] = gene
    return string_protein_gene_map


def get_TF_set():
    from constants import HUMAN_TF_IDLIST_PATH
    tflist = []
    with open(HUMAN_TF_IDLIST_PATH, 'r') as f:
        for line in f.readlines():
            tflist.append(line.strip())
    print(f'processed {len(tflist)} Transcription Factors')
    return tflist 


def construct_HardWiredGenome_A_Matrix(threshold=600, nodes=None):
    """
    Construct the full A matrix of the Hard Wired Genome from the STRING and HURI datasets
    """
    if not os.path.exists(HWG_BASE_PATH):
        os.makedirs(HWG_BASE_PATH)
    print("Constructing Hard Wired Genome A Matrix")

    if nodes is None:
        nodes = get_unique_genes()

    # with open(ENSEMBL_PROTEIN_GENE_PATH, "r") as protein_gene_file:
    #     protein_gene_map = {}
    #     for line in protein_gene_file:
    #         gene, protein = line.split()
    #         protein_gene_map[protein] = gene

    string_protein_gene_map = get_string_protein_gene_map()

    n = len(nodes)

    node_idx = {g: i for i, g in enumerate(nodes)}
    adj = np.zeros((len(nodes), len(nodes)), dtype=int)

    # adj_df = pd.DataFrame(0, index=nodes, columns=nodes)

    # initialize edge set
    edges = set()
    count = 0

    with open(STRING_EXTRACTED_URL_PATH, "r") as f: 
        string_links = f.readlines()

        # skipping header line
        for link in string_links[1:]:
            count += 1
            p1, p2, score = link.strip().split(' ')
            p1 = p1.split('.')[1]
            p2 = p2.split('.')[1]
            score = int(score)
            g1 = string_protein_gene_map.get(p1, None)
            g2 = string_protein_gene_map.get(p2, None)
            # if score >= threshold:
            
            # if g1 in nodes and g2 in nodes:
            if g1 is not None and g2 is not None and score >= threshold:
                edges.add((g1, g2, 1))
    print("STRING interactions processed: ", len(edges))


    with open(HI_UNION_PATH, "r") as huri_links:
        for link in huri_links:
            g1, g2 = link.strip().split()
            edges.add((g1, g2, 1))
    print("HI-UNION interactions processed: ", len(edges))

    with open(LIT_BM_PATH, "r") as lit_links:
        for link in lit_links:
            g1, g2 = link.strip().split()
            edges.add((g1, g2, 1))
    print("LIT-BM interactions processed: ", len(edges))

    tflist = get_TF_set()
    print(tflist[:10])

    # for g1, g2, inter in edges:
    #     if g1 in node_idx and g2 in node_idx:
    #         i, j = node_idx[g1], node_idx[g2]
    #         adj[i, j] = 1
    import pandas as pd

    edgelist = list(edges)
    edgelist_df = pd.DataFrame(edgelist, columns=['source', 'target', 'weight'])

    
    # adj_df = pd.DataFrame(adj, index=nodes, columns=nodes)

    save_path = HARD_WIRED_GENOME_A_EDGELIST.strip('.csv') + f"_{threshold}.csv"
    edgelist_df.to_csv(save_path)
    A_matrix_edgelist = edgelist_df

    print("Constructing Hard Wired Genome B Matrix")
    B_matrix_edgelist = A_matrix_edgelist[
    A_matrix_edgelist['source'].isin(flat_col_order) & A_matrix_edgelist['target'].isin(tflist)
    ]
    print(f'processed edges {len(B_matrix_edgelist)}')

    print("Constructing Hard Wired Genome C Matrix")
    C_matrix_edgelist = A_matrix_edgelist[
        A_matrix_edgelist['source'].isin(tflist) & A_matrix_edgelist['target'].isin(tflist)
    ]
    print(f'processed edges {len(C_matrix_edgelist)}')

    
    save_path = HARD_WIRED_GENOME_B_EDGELIST.strip('.csv') + f"_{threshold}.csv"
    B_matrix_edgelist.to_csv(save_path, index=False)
    
    save_path = HARD_WIRED_GENOME_C_EDGELIST.strip('.csv') + f"_{threshold}.csv"
    C_matrix_edgelist.to_csv(save_path, index=False)

    print(f"Saved Hard Wired Genome to {save_path}")
    return A_matrix_edgelist, B_matrix_edgelist, C_matrix_edgelist

def fetch_HWG(threshold):

    save_path = HARD_WIRED_GENOME_A_EDGELIST.strip('.csv') + f"_{threshold}.csv"
    A_matrix_edgelist = pd.read_csv(save_path)
    
    save_path = HARD_WIRED_GENOME_B_EDGELIST.strip('.csv') + f"_{threshold}.csv"
    B_matrix_edgelist = pd.read_csv(save_path)
    
    save_path = HARD_WIRED_GENOME_C_EDGELIST.strip('.csv') + f"_{threshold}.csv"
    C_matrix_edgelist = pd.read_csv(save_path)

def get_sorted_gene_order():
    # unique_genes = get_unique_genes()
    filtered_gene_metadata = pd.read_csv('./data/filtered_gene_metadata.csv')
    
    nodes = list(filtered_gene_metadata.GeneStableID)
    
    group_dict = filtered_gene_metadata.groupby('ChromosomescaffoldName')['GeneStableID'].apply(list).to_dict()
    
    sorted_groups = dict(sorted(group_dict.items()))
    flat_col_order = [col for group in sorted_groups.values() for col in group]
    
    print(f"ordering {len(flat_col_order)} genes")
    return flat_col_order

def get_TF_lists(log=True):

    adata= anndata.read_h5ad(f'{OPERATIONS_DIRECTORY}/DSET051_B_Matrix_TFsONLY_TFTGDb_TFLink_ChEA3.h5ad')
    
    tfs = adata.obs
    tflist = tfs['TFStableID'].to_list()
    
    repressor = tfs[tfs['TFClass']=='Repressor']
    activator = tfs[tfs['TFClass']=='Activator']
    conflicted = tfs[tfs['TFClass']=='Conflicted']
    
    repressorlist = repressor.TFStableID.to_list()
    activatorlist = activator.TFStableID.to_list()
    conflictedlist = conflicted.TFStableID.to_list()
    if log:
        print(f' Total Activators : {len(activatorlist)}')
        print(f'Total Reprossors : {len(repressorlist)}')
        print(f'Total Conflicted : {len(conflictedlist)}')
    
        print(f'Total Transcription Factors : {len(tflist)}')
    return repressorlist, activatorlist, conflictedlist, tflist

def get_master_regulator_list():
    b_matrix = anndata.read_h5ad(HWG_B_MATRICES)

    bobs = b_matrix.obs_names
    tf_list = b_matrix.obs['TFStableID'].copy()
    
    obs_idx = [b_matrix.obs_names.get_loc(name) for name in bobs]
    var_idx = [b_matrix.var_names.get_loc(name) for name in bobs]
    
    matrix = b_matrix.X
    
    diagonal_values = matrix[obs_idx, var_idx]
    ones_mask = diagonal_values == 1
    
    master_regulator_list = list(b_matrix[bobs[ones_mask]].obs['TFStableID'])
    print(f'total master regulators {len(master_regulator_list)}')
    return master_regulator_list

def construct_HardWiredGenome_B_Matrix(A_Matrix_path=HARD_WIRED_GENOME_A_MATRIX_PATH):

    TF_set = get_TF_set()


    A_matrix = pd.read_csv(A_Matrix_path, index_col=0)

    present_genes = list(set(TF_set) & set(A_matrix.index) & set(A_matrix.columns))

    print(f"Total present genes in A matrix: {len(present_genes)}")
    print(f"Total missing TFs in A matrix: {len(present_genes) - len(TF_set)}")
    B_matrix = A_matrix[present_genes]
    # B_matrix = A_matrix.loc[present_genes, present_genes]

    print(B_matrix.head())
    B_matrix.to_csv(HARD_WIRED_GENOME_B_MATRIX_PATH, index=True)



    # A_matrix = pd.read_csv(A_Matrix_path, index_col=0)

    # transcription_factors = set()

def binarize(mat, thresh):
    if sp.issparse(mat):
        out = mat.copy()
        out.data = (out.data > thresh).astype(np.int8)
        return out
    else:
        return (mat > thresh).astype(np.int8)

def get_a_matrix_threshold(threshold):
    repressorlist, activatorlist, conflictedlist, tf_list = get_TF_lists()
    
    a_matrix_adata = anndata.read_h5ad('/nfs/turbo/umms-indikar/shared/projects/HWG/data/HWG/data/HWG/A_matrix.h5ad')
    a_matrix_adata.layers['HWG_0'] = a_matrix_adata.X
    a_bin = binarize(a_matrix_adata.layers['STRING'], threshold)
    b_bin = a_matrix_adata.layers['HURI']

    merged = a_bin + b_bin 

    if sp.issparse(merged):
        merged.data = np.ones_like(merged.data, dtype=np.int8)
    else:
        merged = (merged > 0).astype(np.int8)

    a_matrix_adata.layers[f'HWG_{threshold}'] = merged
    a_matrix_adata.X = merged

    b_matrix_true = a_matrix_adata[a_matrix_adata.obs_names.isin(tf_list)]

    # c_matrix_true = a_matrix_adata[a_matrix_adata.obs_names.isin(tf_list), a_matrix_adata.var_names.isin(tf_list)]

    c_matrix_true = a_matrix_adata[tf_list, tf_list]
    
    return a_matrix_adata, b_matrix_true, c_matrix_true

def load_htf_motifs():

    filename = TF_MOTIF_PATH
    
    df = pd.read_csv(filename, sep='\t', engine='python', dtype=str)
    
    df.columns = df.columns.str.strip()
    
    df = df.fillna('')  # Or use `NaN` if preferred
    
    return df

from constants import TF_MOTIF_BASE_PATH


def load_cisbp_motifs():

    CISBP_MOTIFS_DIR = f'{TF_MOTIF_BASE_PATH}/Homo_sapiens_2025_07_21_6_25_pm/pwms_all_motifs'
    CISBP_INFO = f'{TF_MOTIF_BASE_PATH}/Homo_sapiens_2025_07_21_6_25_pm/TF_Information_all_motifs.txt'
    print(CISBP_MOTIFS_DIR)


    filename = CISBP_INFO
    
    df = pd.read_csv(filename, sep='\t', engine='python', dtype=str)
    
    df.columns = df.columns.str.strip()
    
    df = df.fillna('')  # Or use `NaN` if preferred

    df_cleaned_first = df.drop_duplicates(subset=['TF_Name'], keep='first')
    
    return df_cleaned_first

    cisbp_motifs = load_cisbp_motifs()


def load_consensus(gene,  scale=1000, gene_id=None, motif_source='cisbp'):

    if gene_id == None:
        gene_id_name_map, gene_name_id_map = generate_gene_id_name_map()
    
        gene_id = gene_name_id_map.get(gene)

    if motif_source == 'htf':
        
        motif_df = load_htf_motifs()
    
        try:
            
            gene_cisbp_id = motif_df[((motif_df['Ensembl ID'] == gene_id) & (motif_df['Best Motif(s)? (Figure 2A)'] == 'TRUE'))]['CIS-BP ID'].values[0]
            gene_motif_path = f'{HTF_MOTIFS_DIR}/{str(gene_cisbp_id)}.txt'
            
        except:
            print('unable to find motif for gene in HTF')
            return None
            
    elif motif_source == 'cisbp':
        motif_df = load_cisbp_motifs()
        gene_cisbp_id = motif_df[motif_df['TF_Name'] == 'TFAP2B']['Motif_ID'].loc[0]
        
        gene_motif_path = f'{HTF_MOTIFS_DIR}/{str(gene_cisbp_id)}.txt'
        
    
    print(f'Gene CIS-BP ID : {gene_cisbp_id}')
    
    gene_motif_path = f'{CISBP_MOTIFS_DIR}/{str(gene_cisbp_id.split("_")[0])}_3.00.txt'
    
    print(gene_motif_path)
    
    gene_pwm = pd.read_csv(gene_motif_path, sep='\t')
    
    # NOTE: THese vales are probabilites so PPM matrices where each column adds up to one
    
    df = gene_pwm
    counts = {
        nt: (df[nt].values * scale).astype(int).tolist()
        for nt in ["A", "C", "G", "T"]
    }
    
    motif = motifs.Motif(counts=counts)
    motif.name = gene
    motif.alphabet = "ACGT"
    motif.pseudocounts = 0.1
    
    consensus = motif.consensus
    print("Consensus sequence:", consensus)
    return consensus, motif, gene_pwm


def setup():
    fetch_ensembl_mapping()
    generate_protein_gene_mappings()
    generate_protein_gene_mapping_STRING()
    generate_unique_gene_set_Huri()
    generate_human_tf_set()
    construct_HardWiredGenome_A_Matrix()
    construct_HardWiredGenome_B_Matrix()


if __name__ == "__main__":
    # setup()
    generate_string_mappings()
    # get_unique_genes()

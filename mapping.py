from biomart import BiomartServer
from constants import DATA_DIRECTORY,  STRING_ALIAS_PATH, STRING_URL_PATH, ENSEMBL_PROTEIN_GENE_PATH, STRING_EXTRACTED_ALIAS_PATH, STRING_EXTRACTED_URL_PATH, ENSEMBL_MAPPING_PATH, UNIQUE_GENE_SET, \
    HURI_URL_PATH, HI_UNION_PATH, LIT_BM_PATH, HUMAN_TF_PATH, HARD_WIRED_GENOME_A_MATRIX_PATH, STRING_PROTEIN_GENE_PATH, STRING_UNIQUE_GENE_SET, HURI_UNIQUE_GENE_SET, HUMAN_TF_SET_PATH, \
    HARD_WIRED_GENOME_B_MATRIX_PATH, HWG_BASE_PATH, ENSEMBL_BASE_PATH

import re
import os
import numpy as np
import pandas as pd

def fetch_ensembl_mapping():
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
            print(f'running line {count}')
            p1, p2, score = line.strip().split(' ')
            lineset = lineset | {p1, p2}

     with open(ENSEMBL_PROTEIN_GENE_PATH, "w") as f:
        for item in lineset:
            f.write(f"{item}\n")
        
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
                unique_genes.add(gene)
                outfile.write(' '.join(matches) + '\n')

    print(f"Total unique genes STRING: {len(unique_genes)}")
    with open(STRING_UNIQUE_GENE_SET, "w") as gene_file:
        for item in unique_genes:
            gene_file.write(f"{item}\n")

        
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

    with open(STRING_UNIQUE_GENE_SET, "r") as f:
        string_nodes = set(line.strip() for line in f)

    with open(UNIQUE_GENE_SET, "r") as f:
        ensemble_nodes = set(line.strip() for line in f)

    with open(HURI_UNIQUE_GENE_SET, "r") as f:
        huri_nodes = set(line.strip() for line in f)

    nodes = string_nodes | ensemble_nodes | huri_nodes
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

def construct_HardWiredGenome_A_Matrix(threshold=600):
    """
    Construct the full A matrix of the Hard Wired Genome from the STRING and HURI datasets
    """
    if not os.path.exists(HWG_BASE_PATH):
        os.makedirs(HWG_BASE_PATH)
    print("Constructing Hard Wired Genome A Matrix")

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

    with open(STRING_EXTRACTED_URL_PATH, "r") as f: 
        string_links = f.readlines()

        # skipping header line
        for link in string_links[1:]:
            p1, p2, score = link.strip().split(' ')
            p1 = p1.split('.')[1]
            p2 = p2.split('.')[1]
            score = int(score)
            g1 = string_protein_gene_map[p1]
            g2 = string_protein_gene_map[p2]
            if score >= threshold:
                edges.add((g1, g2))
    print("STRING interactions processed")


    with open(HI_UNION_PATH, "r") as huri_links:
        for link in huri_links:
            g1, g2 = link.strip().split()
            edges.add((g1, g2))
    print("HI-UNION interactions processed")

    with open(LIT_BM_PATH, "r") as lit_links:
        for link in lit_links:
            g1, g2 = link.strip().split()
            edges.add((g1, g2))
    print("LIT-BM interactions processed")
    
    for g1, g2 in edges:
        if g1 in node_idx and g2 in node_idx:
            i, j = node_idx[g1], node_idx[g2]
            adj[i, j] = 1
    
    adj_df = pd.DataFrame(adj, index=nodes, columns=nodes)

    adj_df.to_csv(HARD_WIRED_GENOME_A_MATRIX_PATH)

    print(f"Saved Hard Wired Genome to {HARD_WIRED_GENOME_A_MATRIX_PATH}")


def construct_HardWiredGenome_B_Matrix(A_Matrix_path=HARD_WIRED_GENOME_A_MATRIX_PATH):

    TF_set = set()
    with open(HUMAN_TF_SET_PATH, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            gene, istf = line.strip().split(',')
            if bool(istf):
                TF_set.add(gene)

    print(f"Total TFs: {len(TF_set)}")


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

if __name__ == "__main__":
    # fetch_ensembl_mapping()
    # generate_protein_gene_mappings()
    # generate_protein_gene_mapping_STRING()
    # generate_unique_gene_set_Huri()
    # generate_human_tf_set()
    # construct_HardWiredGenome_A_Matrix()
    # construct_HardWiredGenome_B_Matrix()

    get_unique_genes()

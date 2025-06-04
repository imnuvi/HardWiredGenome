from biomart import BiomartServer
from constants import DATA_DIRECTORY,  STRING_ALIAS_PATH, STRING_URL_PATH, ENSEMBL_PROTEIN_GENE_PATH, STRING_EXTRACTED_ALIAS_PATH, STRING_EXTRACTED_URL_PATH, ENSEMBL_MAPPING_PATH, UNIQUE_GENE_SET, \
    HURI_URL_PATH, HI_UNION_PATH, LIT_BM_PATH, HUMAN_TF_PATH, HARD_WIRED_GENOME_A_MATRIX_PATH
import re
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

    print(f"Total unique genes: {len(unique_genes)}")
    with open(UNIQUE_GENE_SET, "w") as gene_file:
        for item in unique_genes:
            gene_file.write(f"{item}\n")

def construct_HardWiredGenome_A_Matrix(threshold=600):
    """
    Construct the full A matrix of the Hard Wired Genome from the STRING and HURI datasets
    """
    with open(UNIQUE_GENE_SET, "r") as f:
        nodes = list(line.strip() for line in f)

    with open(ENSEMBL_PROTEIN_GENE_PATH, "r") as protein_gene_file:
        protein_gene_map = {}
        for line in protein_gene_file:
            gene, protein = line.split()
            protein_gene_map[protein] = gene

    n = len(nodes)

    adj_df = pd.DataFrame(0, index=nodes, columns=nodes)

    with open(STRING_EXTRACTED_URL_PATH, "r") as f: 
        string_links = f.readlines()

        # skipping header line
        for link in string_links[1:]:
            p1, p2, score = link.strip().split(' ')
            g1 = protein_gene_map[p1]
            g2 = protein_gene_map[p2]
            if score >= threshold:
                adj_df.at[g1, g2] = 1
    print("STRING interactions processed")

    with open(HI_UNION_PATH, "r") as huri_links:
        for link in huri_links:
            g1, g2 = link.strip().split()
            adj_df.at[g1, g2] = 1

    with open(LIT_BM_PATH, "r") as lit_links:
        for link in lit_links:
            g1, g2 = link.strip().split()
            adj_df.at[g1, g2] = 1

    adj_df.to_csv(HARD_WIRED_GENOME_A_MATRIX_PATH, index=True)

    print(f"Saved Hard Wired Genome to {HARD_WIRED_GENOME_A_MATRIX_PATH}")







if __name__ == "__main__":
#     fetch_ensembl_mapping()
    # generate_string_mappings()
    # generate_protein_gene_mappings()
    construct_HardWiredGenome_A_Matrix()

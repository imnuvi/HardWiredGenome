
DATA_DIRECTORY = '/scratch/indikar_root/indikar1/shared_data/HWG/data'

# STRING data URL
STRING_BASE_PATH = f"{DATA_DIRECTORY}/STRING"

STRING_URL = 'https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz'
STRING_URL_PATH = f"{DATA_DIRECTORY}/STRING/9606.protein.links.v12.0.txt.gz"
STRING_EXTRACTED_URL_PATH = f'{DATA_DIRECTORY}/STRING/9606.protein.links.v12.0.txt'

STRING_ALIAS = 'https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz'
STRING_ALIAS_PATH = f"{DATA_DIRECTORY}/STRING/9606.protein.aliases.v12.0.txt.gz"
STRING_EXTRACTED_ALIAS_PATH = f'{DATA_DIRECTORY}/STRING/9606.protein.aliases.v12.0.txt'

STRING_PROTEIN_LIST_URL = 'https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz'
STRING_PROTEIN_LIST_URL_PATH = f'{DATA_DIRECTORY}/STRING/9606.protein.info.v12.0.txt.gz'
STRING_EXTRACTED_PROTEIN_LIST_URL_PATH = f'{DATA_DIRECTORY}/STRING/9606.protein.info.v12.0.txt'


STRING_PROTEIN_GENE_PATH = f"{DATA_DIRECTORY}/STRING/protein_gene_map.txt"
STRING_UNIQUE_GENE_SET = f"{DATA_DIRECTORY}/STRING/unique_gene_set.txt"


# HURI data URLs
HURI_BASE_PATH = f"{DATA_DIRECTORY}/HURI"
HURI_URL = 'http://www.interactome-atlas.org/data/HuRI.tsv'
HURI_URL_PATH = f"{DATA_DIRECTORY}/HURI/HuRI.tsv"

HI_UNION = 'http://www.interactome-atlas.org/data/HI-union.tsv'
HI_UNION_PATH = f"{DATA_DIRECTORY}/HURI/HI-union.tsv"

LIT_BM = 'http://www.interactome-atlas.org/data/Lit-BM.tsv'
LIT_BM_PATH = f"{DATA_DIRECTORY}/HURI/Lit-BM.tsv"

HURI_UNIQUE_GENE_SET = f"{DATA_DIRECTORY}/HURI/unique_gene_set.txt"

# HUMAN TF data URL
HUMAN_TF_BASE_PATH = f"{DATA_DIRECTORY}/HTF"
HUMAN_TF = 'https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv'
HUMAN_TF_PATH = f"{DATA_DIRECTORY}/HTF/DatabaseExtract_v_1.01.csv"

HUMAN_TF_SET_PATH = f"{DATA_DIRECTORY}/HTF/human_tf_set.txt"


# Ensembl mapping data
ENSEMBL_BASE_PATH = f"{DATA_DIRECTORY}/ENSEMBL"
ENSEMBL_MAPPING_PATH = f"{DATA_DIRECTORY}/ENSEMBL/ensembl_mapping.txt"
ENSEMBL_PROTEIN_GENE_PATH = f"{DATA_DIRECTORY}/ENSEMBL/protein_gene_map.txt"
UNIQUE_GENE_SET = f"{DATA_DIRECTORY}/ENSEMBL/unique_gene_set.txt"


HWG_BASE_PATH = f"{DATA_DIRECTORY}/HWG"
HARD_WIRED_GENOME_A_MATRIX_PATH = f"{DATA_DIRECTORY}/HWG/A_Matrix.csv"
HARD_WIRED_GENOME_B_MATRIX_PATH = f"{DATA_DIRECTORY}/HWG/B_Matrix.csv"
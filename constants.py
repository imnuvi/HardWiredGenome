
# BASEDIR = '/scratch/indikar_root/indikar1/shared_data/HWG'

BASEDIR = '/nfs/turbo/umms-indikar/shared/projects/HWG/data/HWG'
# BASEDIR = '/Users/ramprakash/development/lab_projects/Rajapakse_lab/data/HWG'

DATA_DIRECTORY = f'{BASEDIR}/data'
OPERATIONS_DIRECTORY = f'{BASEDIR}/operations'

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


STRING_PROTEIN_SET_FROM_LINKS = f"{DATA_DIRECTORY}/STRING/protein_set_from_links.txt"
STRING_PROTEIN_GENE_PATH = f"{DATA_DIRECTORY}/STRING/protein_gene_map.txt"
STRING_UNIQUE_PROTEIN_SET = f"{DATA_DIRECTORY}/STRING/unique_protein_set.txt"
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


HUMAN_TF_IDLIST_URL = 'https://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt'
HUMAN_TF_IDLIST_PATH = f"{DATA_DIRECTORY}/HTF/TFs_Ensembl_v_1.01.txt"


HUMAN_TF_SET_PATH = f"{DATA_DIRECTORY}/HTF/human_tf_set.txt"


# Ensembl mapping data
ENSEMBL_BASE_PATH = f"{DATA_DIRECTORY}/ENSEMBL"
ENSEMBL_MAPPING_PATH = f"{DATA_DIRECTORY}/ENSEMBL/ensembl_mapping.txt"
ENSEMBL_PROTEIN_GENE_PATH = f"{DATA_DIRECTORY}/ENSEMBL/protein_gene_map.txt"
UNIQUE_GENE_SET = f"{DATA_DIRECTORY}/ENSEMBL/unique_gene_set.txt"


HWG_BASE_PATH = f"{DATA_DIRECTORY}/HWG"
HARD_WIRED_GENOME_A_EDGELIST = f"{DATA_DIRECTORY}/HWG/A_Matrix_edgelist.csv"
HARD_WIRED_GENOME_B_EDGELIST = f"{DATA_DIRECTORY}/HWG/B_Matrix_edgelist.csv"
HARD_WIRED_GENOME_C_EDGELIST = f"{DATA_DIRECTORY}/HWG/C_Matrix_edgelist.csv"
HARD_WIRED_GENOME_A_MATRIX_PATH = f"{DATA_DIRECTORY}/HWG/A_Matrix.csv"
HARD_WIRED_GENOME_B_MATRIX_PATH = f"{DATA_DIRECTORY}/HWG/B_Matrix.csv"
HARD_WIRED_GENOME_C_MATRIX_PATH = f"{DATA_DIRECTORY}/HWG/C_Matrix.csv"
HARD_WIRED_GENOME_BP_MATRIX = f"{DATA_DIRECTORY}/HWG/BP_Matrix.csv"


HWG_A_MATRICES = f"{HWG_BASE_PATH}/A_matrix.h5ad"
HWG_B_MATRICES = f'{OPERATIONS_DIRECTORY}/B_matrices.h5ad'

RNASEQ_PATH = f'{DATA_DIRECTORY}/RNAseq'

OPERATIONS_PATH = f'{DATA_DIRECTORY}/operations'

TF_MOTIF_BASE_PATH = f'{DATA_DIRECTORY}/TF_MOTIFS'
TF_MOTIF_PATH = f'{DATA_DIRECTORY}/TF_MOTIFS/Human_TF_MotifList_v_1.01.txt'

HTF_MOTIFS_URL = 'https://humantfs.ccbr.utoronto.ca/download/v_1.01/PWMs.zip'
HTF_MOTIFS_PATH = f'{TF_MOTIF_BASE_PATH}/PWMs.zip'
HTF_MOTIFS_DIR = f'{TF_MOTIF_BASE_PATH}/HTF_PWMs'  # Destination folder


# Reference genomic data URLs
REFERENCE_GTF_HG38_URL = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
REFERENCE_DIR = f'{DATA_DIRECTORY}/REFERENCE'
REFERENCE_GTF_HG38_PATH = f'{REFERENCE_DIR}/gencode.v46.annotation.gtf.gz.feather'

REFERENCE_GENES_BED_URL = 'https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2773508654_C64YruQsla02ZUjIh3e4Q1gqq5rq&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED'
REFERENCE_GENES_BED_PATH = f'{REFERENCE_DIR}/genes.bed'

REFERENCE_CHROM_SIZES_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
REFERENCE_CHROM_SIZES_PATH = f'{REFERENCE_DIR}/hg38.chrom.sizes'

REFERENCE_GENOME_URL = 'https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz'
REFERENCE_GENOME_PATH = f'{REFERENCE_DIR}/hg38.fa'

GENCODE_ANNOTATION_URL = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz'
GENCODE_ANNOTATION_PATH = f'{REFERENCE_DIR}/gencode.v43.annotation.gtf.gz'


# NOTE: CIS-BP has to be downloaded manually from here - https://cisbp.ccbr.utoronto.ca/bulk.php. Unzip and add a constant here if needed.

def modify_HWG_path(suffix):
    save_path = HARD_WIRED_GENOME_A_MATRIX_PATH.strip('.csv') + f"_{suffix}.csv"
    return save_path


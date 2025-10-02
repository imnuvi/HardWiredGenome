
import seqlogo
import numpy as np
import pandas as pd
import anndata


import subprocess
import os

from Bio import motifs
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.motifs import meme


from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.data import ontology
from alphagenome.visualization import plot_components
from alphagenome.models import dna_client



from mapping import get_sorted_gene_order, generate_gene_id_name_map, get_master_regulator_list, get_TF_lists, get_a_matrix_threshold, load_htf_motifs, generate_gene_id_name_map, load_consensus

from genome_utils import write_meme_format, parse_fasta, write_df_to_meme

from constants import API_KEY, HTF_MOTIFS_DIR, CISBP_MOTIFS_DIR, REFERENCE_GENOME_PATH, GENCODE_ANNOTATION_PATH, REFERENCE_GTF_HG38_PATH, REFERENCE_DIR, TF_MOTIF_BASE_PATH, REFERENCE_GENOME_PATH

master_regulator_list = get_master_regulator_list()
repressorlist, activatorlist, conflictedlist, tf_list = get_TF_lists()

gene_id_name_map, gene_name_id_map = generate_gene_id_name_map()



tempdir = '/nfs/turbo/umms-indikar/shared/projects/HWG/data/HWG/temp'

b_matrix = anndata.read_h5ad('/nfs/turbo/umms-indikar/shared/projects/HWG/data/HWG/operations/B_matrices_true.h5ad')

genes = b_matrix.var_names
TFs = b_matrix.obs_names


gene_counts = len(genes)
TF_counts = len(TFs)


print('TFs : ', gene_counts)
print('Genes : ' , TF_counts)


# helper functions

        
import matplotlib.pyplot as plt
    


def parse_unique_ids(value):
    return value.split('.')[0]

def parse_id_attr(row, keys, operation = parse_unique_ids):
    resdict = {}
    for key in keys:
        resdict[key] = operation(row[key])
        
    return pd.Series(resdict)

def pick_best_transcript(group):
    under_100kb = group[group["length"] <= 100000]
    if not under_100kb.empty:
        return under_100kb.sort_values("length", ascending=False).iloc[0]
    else:
        return group.sort_values("length", ascending=True).iloc[0]
        # group["dist_to_100kb"] = (group["length"] - 100_000).abs()
        # return group.sort_values("dist_to_100kb").iloc[0]



def plot_distribution(read_type, df, field):
    # Checking average gene length 
    
    
    mean_val = df[field].mean()
    median_val = df[field].median()
    max_val = df[field].max()
    
    print(f"Total {read_type}s: {len(df)}")
    print(f"Median {read_type} length: {median_val:,.0f} bp")
    print(f"Mean {read_type} length: {mean_val:,.0f} bp")
    print(f"Max {read_type} length: {max_val:,.0f} bp")
    
    
    plt.figure(figsize=(10,6))
    plt.hist(df[field], bins=100, color='skyblue', edgecolor='black', log=True)
    
    plt.axvline(mean_val, color='r', linestyle='--', label=f'Mean: {mean_val:.2f}')
    plt.axvline(median_val, color='r', linestyle=':', label=f'Median: {median_val:.2f}')

    
    plt.legend()

    plt.xlabel(f'{read_type} length (bp)')
    plt.ylabel('Frequency (log scale)')
    plt.title(f'Distribution of {read_type} Lengths (GENCODE v46)')
    plt.xlim(0, 2000000)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()
    # plt.savefig('figures/{read_type}_{field}s.png')


import pandas as pd
from Bio import SeqIO, motifs
from Bio.Seq import Seq
from pathlib import Path

import math

def fetch_reference_locations():

    gtf = pd.read_feather(REFERENCE_GTF_HG38_PATH)


    
    genes_df = gtf[gtf['Feature'] == 'gene'].copy()
    keys = ['gene_id']
    genes_df[keys] = genes_df.apply(
        parse_id_attr,
        axis=1,
        args=(keys, parse_unique_ids)
    )
    genes_df = genes_df[genes_df['gene_id'].isin(genes)].copy()
    genes_df['length'] = genes_df['End'] - genes_df['Start'] + 1

    
    transcripts = gtf[gtf['Feature'] == 'transcript'].copy()
    keys = ['gene_id', 'gene_type', 'transcript_id']
    transcripts[keys] = transcripts.apply(
        parse_id_attr,
        axis=1,
        args=(keys, parse_unique_ids)
    )
    transcripts = transcripts[transcripts['gene_id'].isin(genes)].copy()
    transcripts['length'] = transcripts['End'] - transcripts['Start'] + 1


    best_transcripts = transcripts.groupby("gene_id", group_keys=False).apply(pick_best_transcript, include_groups=True)

    return genes_df, transcripts, best_transcripts

def get_tss_flank(row, flank=1000):
    if row["Strand"] == "+":
        tss = row["Start"]
        flank_start = max(0, tss - flank)
        flank_end = tss + flank
    else:
        tss = row["End"]
        flank_start = max(0, tss - flank)
        flank_end = tss + flank
    return pd.Series([flank_start, flank_end])

def fetch_bed_files(flank_size=10000):
    # Get the TSS FLanking region

    genes, transcripts, best_transcripts = fetch_reference_locations()

    bed_csv = f"{REFERENCE_DIR}/protein_coding_tss_flank_{flank_size}.bed"
    one_mb_flanks = f"{REFERENCE_DIR}/protein_coding_tss_flank_1mb.bed"
    two_mb_flanks = f"{REFERENCE_DIR}/protein_coding_tss_flank_2mb.bed"




    if not os.path.exists(bed_csv):
        best_transcripts[["flank_start", "flank_end"]] = best_transcripts.apply(get_tss_flank, axis=1, args=(flank_size,))
        bed = best_transcripts[["Chromosome", "flank_start", "flank_end", "gene_id", "Strand", "Score"]].copy()
        bed = bed[["Chromosome", "flank_start", "flank_end", "gene_id", "Score", "Strand"]]
        bed.to_csv(bed_csv, sep="\t", header=False, index=False)
    else:
        bed = pd.read_csv(bed_csv, sep="\t")



    if not os.path.exists(one_mb_flanks):
        best_transcripts[["flank_start", "flank_end"]] = best_transcripts.apply(get_tss_flank, axis=1, args=(500000,))
        
        
        
        bed_one_m = best_transcripts[["Chromosome", "flank_start", "flank_end", "gene_id", "Strand", "Score"]].copy()
        bed_one_m = bed_one_m[["Chromosome", "flank_start", "flank_end", "gene_id", "Score", "Strand"]]
        
        bed_one_m.to_csv(one_mb_flanks, sep="\t", header=False, index=False)
    else:
        bed_one_m= pd.read_csv(one_mb_flanks, sep="\t", header=None)



    if not os.path.exists(two_mb_flanks):
        best_transcripts[["flank_start", "flank_end"]] = best_transcripts.apply(get_tss_flank, axis=1, args=(1000000,))
        
        
        
        bed_two_m = best_transcripts[["Chromosome", "flank_start", "flank_end", "gene_id", "Strand", "Score"]].copy()
        bed_two_m = bed_two_m[["Chromosome", "flank_start", "flank_end", "gene_id", "Score", "Strand"]]
        
        bed_two_m.to_csv(two_mb_flanks, sep="\t", header=False, index=False)
    else:
        bed_two_m= pd.read_csv(two_mb_flanks, sep="\t", header=None)

    return bed_one_m, bed_two_m, bed_csv


def call_alphagenome(sequence, chromosome, interval_start, interval_end, seqtype='generic'):
    
    dna_model = dna_client.create(API_KEY)
    
    
    # output = dna_model.predict_sequence(
    #     sequence='GATTACA'.center(2048, 'N'),  # Pad to valid sequence length.
    #     requested_outputs=[dna_client.OutputType.RNA_SEQ],
    #     ontology_terms=['UBERON:0002048'],  # Lung.
    # )
    
    # interval = genome.Interval(chromosome, interval_start, interval_start + 2048, '+')
    interval = genome.Interval(chromosome, interval_start, interval_end, '+')

    ontologies = ['UBERON:0001114', 'UBERON:0002107', 'UBERON:0002048', 'UBERON:0002113', 'UBERON:0002097']
    # right liver lobe, liver, lung, kidney,  skin
    
    output = dna_model.predict_sequence(
        interval=interval,
        # sequence='GATTATATATATTTTACCAACTTTGGGGGGGATTATTTCCGGGGGGGGCACACCCAAA'.center(2048, 'N'),  # Pad to valid sequence length.
        sequence = sequence,
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
        ontology_terms=ontologies,
    )  # Right liver lobe.
    
    print(output.rna_seq.values.shape)

    

    plot_components.plot(
        components=[
            # plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(output.rna_seq),
        ],
        interval=interval,
    )
    
    # plt.show()
    plt.savefig(f'figures/B_prime/{seqtype}.png', dpi=300, bbox_inches = 'tight')

def compare_output_expression():
    pass
    

def replace_motif(raw_sequence, promoter_region, fimo_matches, original_motif, replacement_motif, seqlen=1000000, method='group', replacement_counts=None):
    """
    Keeps the TSS as the centre of the window and performs replacement along the sides. Each replacement can have, length differences. 
    inputs:
    - raw 1MB sequence region surrounding the TSS
    - 20kb promoter region surrounding the TSS
    - Fimo matches for the background TF
    - Motif for the background TF
    - Motif for the target TF
    """

    raw_len = len(raw_sequence.seq)
    raw_seq = raw_sequence.seq
    raw_chrom, raw_range = raw_sequence.id.split(':')
    raw_start, raw_end  = list(map(int, raw_range.split('-')))

    promoter_seq = promoter_region.seq
    promoter_len = len(promoter_region.seq)
    prom_chrom, prom_range = promoter_region.id.split(':')
    prom_start, prom_end  = list(map(int, prom_range.split('-')))
    
    chromosome = prom_chrom

    print(f'Region Start : {raw_start} and Region End : {raw_end}')
    tss_index = (raw_end - raw_start )//2
    print('TSS of raw seq', tss_index)
    
    tss_index_promoter = (prom_end - prom_start )//2
    print('TSS of promoter seq', tss_index_promoter)

    fimo_matches.sort_values(by='start')

    if replacement_counts:
        fimo_matches = fimo_matches[:1]
        
    new_prom = ''
    temp_pointer = 0

    ordered_matches = []

    print(fimo_matches)
    print('Length of fimo matches', len(fimo_matches))

    
    for index, match in fimo_matches.iterrows():
        if math.isnan(match.start) or math.isnan(match.stop):
            continue
        # print(f'Match Start: {match.start} and Match Stop: {match.stop}')
        index_start = int(match.start - prom_start)
        index_end = int(match.stop - prom_start)

        ### CURRENT ALGORITHM
        ## in case of overlap, replace the first occurance till the next occurance and then replace only the non overlaps, without caring for the other overlaps.
        # This part is plug and play. Experimenting with this.

        if method == 'group':
            if index_start > temp_pointer:
                new_prom += promoter_seq[temp_pointer:index_start]
                temp_pointer = index_end
                new_prom += replacement_motif
            else:
                continue
        elif method == 'full':
            new_prom += replacement_motif

        # print(match.motif_id)

        # To check the replacements if needed

    
        # print(new_prom[index_start-5:index_start], new_prom[index_start:index_end], new_prom[index_end:index_end+5])
        # print(promoter_seq[index_start-5:index_start], promoter_seq[index_start:index_end], promoter_seq[index_end:index_end+5])

    # replacement happens with the promoter start as the center

    prom_start_full = prom_start-raw_start
    prom_end_full = prom_end-raw_start
    base_seq = raw_seq[prom_start_full-(seqlen//2): prom_start_full+(seqlen//2)]

    repl_start_index = prom_start_full + len(new_prom)
    repl_end_index = repl_start_index + (seqlen//2) - len(new_prom)
    print(repl_start_index, repl_end_index, repl_start_index + repl_end_index)

    print(len(new_prom))
    print('---------------------')
    
    replaced_seq = raw_seq[prom_start_full-(seqlen//2):prom_start_full] + new_prom + raw_seq[repl_start_index:repl_end_index]

    # print('*****************************')
    # print(str(base_seq[:100]))
    # print(replaced_seq[:100])

    print('%%%%%%%%%%%%%%')
    print(len(base_seq))
    print(len(replaced_seq))
    
    print('%%%%%%%%%%%%%%')

    # print(prom_start, raw_start)

    
    # print(len(new_prom), len(promoter_seq))
    # print('#$$$$$$$$$$$$$$$$$$$$$$$$$$')
    # print(original_motif)
    # print(replacement_motif)
    
    return base_seq, replaced_seq, prom_start_full-(seqlen//2), prom_start_full+(seqlen//2), chromosome

def replacement_runner(base_gene_id, base_tf_id, target_tf_id):

    base_tf_name = gene_id_name_map.get(base_tf_id)
    base_tf = base_tf_id
    target_name = gene_id_name_map.get(target_tf_id)
    target_id = target_tf_id

    print(f"Running replacement for base TF: {base_tf_name} with target TF: {target_name}")
    print(f"Running replacement for base TF ID: {base_tf} with target TF ID: {target_id}")
    gene_name = gene_id_name_map.get(base_gene_id)
    RAW_SEQ_FASTA = f'{REFERENCE_DIR}/{gene_name}_2mb_region.fa'
    PROMOTER_REGION_FASTA = f'{REFERENCE_DIR}/{gene_name}_promoter.fa'
    
    FIMO_MATCHES = 'fimo_out/fimo.tsv'

    
    
    orig_consensus, orig_motif, orig_pwm = load_consensus(base_tf_name, gene_id=base_tf)
    repl_consensus, repl_motif, repl_pwm = load_consensus(target_name, gene_id=target_id)

    
    fimo_tsv_out = pd.read_csv(FIMO_MATCHES, sep='\t')


    for record in SeqIO.parse(RAW_SEQ_FASTA, "fasta"):
        raw_sequence = record
        
        
    for record in SeqIO.parse(PROMOTER_REGION_FASTA, "fasta"):
        promoter_region = record

    # For alphagenome
    # [2048, 16384, 131072, 524288, 1048576]
    # seqlen = 524288
    seqlen = 16384
    method = 'full'

    base_seq, repl_seq, interval_start, interval_end, chromosome = replace_motif(raw_sequence, promoter_region, fimo_tsv_out, orig_consensus, repl_consensus, seqlen, method)
    base_seq = str(base_seq).upper()

    repl_seq = str(repl_seq).upper()

    

    # start = genes_df[genes_df['gene_name'] == gene_name].iloc[0].loc['Start']
    # end = genes_df[genes_df['gene_name'] == gene_name].iloc[0].loc['End']

    # call_alphagenome(base_seq, chromosome, start, end, f'{gene_name}_{target_name}_{method}_raw')
    # call_alphagenome(repl_seq, chromosome, start, end, f'{gene_name}_{target_name}_{method}_repl')
    

    call_alphagenome(base_seq, chromosome, interval_start, interval_end, f'{gene_name}_{target_name}_{method}_raw')
    call_alphagenome(repl_seq, chromosome, interval_start, interval_end, f'{gene_name}_{target_name}_{method}_repl')
    

# for target_tf_id in valid_repressors[:1]:

# for target_tf_id in [gene_name_id_map.get(name) for name in ['FLI1', 'TAL1']]:
#     replacement_runner(base_gene_id, target_tf_id)


def run_pipeline(gene):

    genedir = f'{tempdir}/{gene}'
    if not os.path.exists(genedir):
        os.makedirs(genedir)

    bed_one_m, bed_two_m, bed_flanks = fetch_bed_files()

    gene_id = gene_name_id_map.get(gene)

    gene_filtered_bed_one_m = bed_one_m[bed_one_m[3] == gene_id]
    gene_filtered_bed_two_m = bed_two_m[bed_two_m[3] == gene_id]
    gene_filtered_flank = bed_flanks[bed_flanks[3] == gene_id]

    gene_one_m = f'{genedir}/{gene}_1mb_flank.bed'
    gene_two_m = f'{genedir}/{gene}_2mb_flank.bed'
    gene_filtered_flank_m = f'{genedir}/{gene}_10kb.bed'

    gene_one_m_fa = f'{genedir}/{gene}_1mb_flank.fa'
    gene_two_m_fa = f'{genedir}/{gene}_2mb_flank.fa'
    gene_filtered_flank_m_fa = f'{genedir}/{gene}_10kb.fa'
    

    gene_filtered_bed_one_m.to_csv(gene_one_m, sep='\t', header=False, index=False)
    gene_filtered_bed_two_m.to_csv(gene_two_m, sep='\t', header=False, index=False)

    # gene_filtered_bed_two_m

    genes_df, transcripts, best_transcripts = fetch_reference_locations()
    print(genes_df, transcripts, best_transcripts)
    print(bed_one_m, bed_two_m)

    ### Running Bedtools to extract the fasta file for the target gene
    print("Running Bedtools and writing output to : ", tempdir)

    try:

        # bedtools getfasta -fi REFERENCE/hg38.fa -bed U2AF2_1mb.bed > U2AF2_1mb_region.fa
        subcommand = f'bedtools getfasta -fi {REFERENCE_GENOME_PATH} -bed {gene_one_m} > {gene_one_m_fa}' 
        result = subprocess.run(subcommand, shell=True, capture_output=True, text=True, check=True)

        subcommand = f'bedtools getfasta -fi {REFERENCE_GENOME_PATH} -bed {gene_two_m} > {gene_two_m_fa}'
        result = subprocess.run(subcommand, shell=True, capture_output=True, text=True, check=True)

        print('Done')

    except Exception as e:
        print("error occured processing bedtools: ", e)


    from network_utils import GeneUtils


    GUtils = GeneUtils()
    gene_influencers = GUtils.get_influencers(gene)

    print(gene_influencers)

    # consensus, motif, gene_pwm = load_consensus(base_tf_name, gene_id=base_tf)

    # savepath = f'{tempdir}/{base_tf}.meme'

    gene_id = 'ENSG00000063244'

    repressorlist, activatorlist, conflictedlist, tf_list = get_TF_lists()
    
    # base_tf = ''
    # target_tf = ''
    # gene_name_id_map.get('')


    # replacement_runner(gene_id, base_tf_id, target_tf_id)



# gene_target = 'ACADM'
# gene_target_id = gene_name_id_map.get(gene_target)
# run_pipeline(gene_target_id)


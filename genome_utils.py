
from Bio import SeqIO


def write_meme_format(motif, filepath, pseudocounts=0.01):
    """
    Write a Biopython motif object to MEME format manually.
    Assumes motif.counts is available (i.e., PFM, not PWM).
    """
    with open(filepath, "w") as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies:\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        
        f.write(f"MOTIF {motif.name}\n")
        f.write("letter-probability matrix: alength= 4 w= {} nsites= 20 E= 0\n".format(motif.length))
        
        for i in range(motif.length):
            line = []
            for base in "ACGT":
                val = motif.counts[base][i] + pseudocounts
                line.append(f"{val:.6f}")
            f.write(" ".join(line) + "\n")

def write_df_to_meme(df, filepath, motif_name="Motif1"):
    """
    Write a MEME motif file from a pandas DataFrame with probabilities.
    
    Parameters:
        df: DataFrame with columns ['Pos', 'A', 'C', 'G', 'T']
        output_file: File path to write the MEME motif
        motif_name: Motif identifier
    """
    # Drop 'Pos' column if present
    if 'Pos' in df.columns:
        df = df.drop(columns='Pos')
    
    # Ensure correct column order
    df = df[['A', 'C', 'G', 'T']]

    # Transpose so rows are bases, columns are positions
    ppm_df = df.T
    width = ppm_df.shape[1]
    
    print(f'writing meme format to {filepath}')

    with open(filepath, 'w') as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies (from uniform background):\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")

        f.write(f"MOTIF {motif_name}\n")
        f.write(f"letter-probability matrix: alength= 4 w= {width} nsites= 20 E= 0\n")

        for col in ppm_df.columns:
            row = [f"{ppm_df.loc[base, col]:.6f}" for base in ['A', 'C', 'G', 'T']]
            f.write(" ".join(row) + "\n")

    

def parse_fasta(fasta_path, limit=None):

    it = 0
    recordlist = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        if limit and it == limit:
            break
        recordlist.append(record)
        print(record.id)
        print(len(record.seq))
        it += 1
    
    return recordlist
    
import os
import time
import subprocess
import glob
import shutil
import pandas as pd
import numpy as np
import shlex
from Bio import SeqIO

from pipeline.useful_funct import contig_to_genome
from pipeline.metadata_info import add_genome_metadata

# --------Definition of functions
#------------------------ Function to run blastn
def run_blastn(query, output, db, max_target, evalue_threshold = 10):
    """
    Executes the blastn command with the specified query, output file, and database.

    Args
    ------
    query (str): The path to the query fasta file.
    output (str): The path to the output file where results will be saved.
    db (str): The path to the blast database.
    max_target (int): Maximum number of aligned sequences keeped into the results.
    evalue_threshold (int): E-value threshold.
                            Default to 10.

    Return
    ------
    output (str): The path to the output file where results are saved.
    elapsed (float): The time taken to execute the blastn command in seconds.
    """
    # Creates the output directory to avoid errors
    os.makedirs(os.path.dirname(output), exist_ok=True)
    # Blastn command
    blast_cmd = ["blastn", "-query", query, # Query fasta file
                 "-db", db, # Blast database
                 "-out", output, 
                 "-max_target_seqs", str(max_target), # Max number of aligned seqs to keep
                 "-evalue", str(evalue_threshold),
                 "-outfmt", "6"] # Output format 6 (tabular)
    start = time.time() # Start time to measure execution time
    # Execute the command and handle errors
    print(f"Executing: \n{shlex.join(blast_cmd)}")
    try:
        subprocess.run(blast_cmd, 
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing blastn:\n{e.stderr}")
        return

    total_time = time.time() - start # End time to measure execution time

    return output, total_time

def process_blast_tsv(path):
    """
    Reads a BLAST tabular output (outfmt 6) and adds strand and normalized
    subject coordinates.

    Args
    ------
    path (str): The path to the BLAST tabular output file.

    Return
    ------
    df (pd.DataFrame): One row per hit, with the 12 standard outfmt 6 columns plus:
                       - strand (str): "+" if seq_start <= seq_end, "-" otherwise.
                       - start (int): Lowest subject coordinate of the hit.
                       - end (int): Highest subject coordinate of the hit.
                       
                       If the file is empty, an empty DataFrame with the same
                       columns is returned.
    """
    columns = ["query_id", "seq_id", "identity", "length", "mismatch", "gapopen", 
               "query_start", "query_end", "seq_start", "seq_end", "Evalue", "BitScore"]
    try: # Reads the csv file
        df = pd.read_csv(path, sep="\t", names = columns)
    except pd.errors.EmptyDataError:
        # If no hits,, keeps the expected columns 
        df = pd.DataFrame(columns = columns)
    # If start is lower than end, the strand is positive
    pos_strand = df["seq_start"] <= df["seq_end"]
    # Assign the strand column based on the comparison of seq_start and seq_end
    df["strand"] = np.where(pos_strand, "+", "-")  
    # Assign the start and end columns based on the minimum and maximum of seq_start and seq_end
    df["start"] = df[["seq_start", "seq_end"]].min(axis=1)
    df["end"] = df[["seq_start", "seq_end"]].max(axis=1)
    return df

def extract_hit_sequences(hits_df):
    """
    Extracts the nucleotide sequence of each BLAST hit from its genome fasta file.
    
    Args
    ------
    hits_df (pd.DataFrame): Hits table with, at least, the columns
                            genome_dir, seq_id, start, end and strand.

    Return
    ------
    hits_df (pd.DataFrame): A copy of the input with an additional "sequence" column. 
    """
    hits_df = hits_df.copy()
    hits_df["sequence"] = ""

    # Group by genome so each .fna file is read only once
    for genome_dir, group in hits_df.groupby("genome_dir"):
        # Get the path to the .fna file in the genome directory
        fna_ext = glob.glob(os.path.join(genome_dir, "*.fna"))
        fna = fna_ext[0]

        # Parse the fasta file and create a dictionary of contig sequences
        contigs = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))

        for idx, row in group.iterrows(): # Iterate over each hit in the group
            # Extract the sequence of the hit from the contig sequence
            # (BLAST coordinates are 1-based and inclusive)
            seq = contigs[row["seq_id"]].seq[row["start"] - 1 : row["end"]]
            if row["strand"] == "-":
                seq = seq.reverse_complement()  # Orient the sequence like the query
            # Add the extracted sequence to the hits_df
            hits_df.at[idx, "sequence"] = str(seq)

    return hits_df

def process_hits(tsv_path, dir_genomes, multifasta, extra_columns=None):
    """
    Parses the BLAST results, removes the genomes without hits and extracts the
    sequences of the remaining hits.

    Args
    ------
    tsv_path (str): The path to the BLAST tabular output file.
    dir_genomes (str): The path to the directory containing one subdirectory per genome.
    multifasta (str): The path to the multifasta file used as BLAST database source.
    extra_columns (dict, optional): Extra columns to add, as
                                    {column name: dotted path in the JSON}.
                                    Defaults to None (only species is added).

    Return
    ------
    final_df (pd.DataFrame): One row per hit with the columns genome, genome_dir, seq_id, 
                             start, end, strand and sequence.

    """
    # Read the BLAST results and add strand and normalized coordinates
    hits_df = process_blast_tsv(tsv_path)  
    contig2genome = contig_to_genome(dir_genomes)

    # Control to prevent deleting everything by mistake
    if hits_df.empty: # If BLAST returned no hits, something failed
        raise ValueError("BLAST returned no hits: nothing was deleted.")
    # If there are seq_ids without a corresponding genome, something failed
    if not hits_df["seq_id"].isin(contig2genome).all(): 
        raise ValueError("Hay seq_id sin genoma asignado: no borro nada.")
    
    # Add genome_dir and genome columns to hits_df
    hits_df["genome_dir"] = hits_df["seq_id"].map(contig2genome)
    hits_df["genome"] = hits_df["genome_dir"].map(os.path.basename)

    # Genomes without any hit are removed from disk
    genomes_with_hit = set(hits_df["genome_dir"])
    to_delete = sorted(set(contig2genome.values()) - genomes_with_hit)
    for each_del in to_delete:
        shutil.rmtree(each_del)
    # Removes multifasta file
    os.remove(multifasta)
    
    df = hits_df[["genome", "genome_dir", "seq_id", "start", "end", "strand"]]
    # Extracts the nucleotide sequence of each hit from its genome fasta file
    final_df = extract_hit_sequences(df)
    final_df = add_genome_metadata(final_df, extra_columns)

    return final_df

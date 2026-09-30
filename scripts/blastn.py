import os
import subprocess
import glob
import shutil
import pandas as pd
import numpy as np
import shlex
from Bio import SeqIO

# --------Definition of functions
#------------------------ Function to run blastn
def run_blastn(query, output, db, max_target, evalue_threshold = 10, threads = 18):
    """
    Executes the blastn command with the specified query, output file, and database.

    Args
    ------
    query (str): The path to the query fasta file.
    output (str): The path to the output file where results will be saved.
    db (str): The path to the blast database.
    max_target (int): Maximum number of aligned sequences keeped into the results.
    evalue_threshold (int, optional): E-value threshold.
                                      Default to 10.
    threads (int, optional): Number of threads to use.
                             Default to 18.

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
                 "-num_threads", str(threads),
                 "-outfmt", "6"] # Output format 6 (tabular)
    # Execute the command and handle errors
    print(f"Executing: \n{shlex.join(blast_cmd)}")
    try:
        subprocess.run(blast_cmd, 
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing blastn:\n{e.stderr}")
        raise TimeoutError # to avoid removing all the genomes

    return output

def process_blast_tsv(path):
    columns = ["query_id", "seq_id", "identity", "length",
           "mismatch", "gapopen", "query_start", "query_end",
           "seq_start", "seq_end", "Evalue", "BitScore"]
    try:
        df = pd.read_csv(path, sep="\t", names = columns)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns = columns)

    df["strand"] = np.where(df["seq_start"] <= df["seq_end"], "+", "-")
    df["start"] = df[["seq_start", "seq_end"]].min(axis=1)
    df["end"] = df[["seq_start", "seq_end"]].max(axis=1)
    return df

def contig_to_genome(dir_genomes):
    contig2genome = {}
    pattern = os.path.join(dir_genomes, "*", "*", "ncbi_dataset", "data", "*", "*.fna")
    for fna in glob.glob(pattern):
        genome_path = os.path.dirname(fna)  # .../data/GCF_xxxxx
        with open(fna) as f:
            for line in f:
                if line.startswith(">"):
                    contig2genome[line[1:].split()[0]] = genome_path
    return contig2genome



def extract_hit_sequences(hits_df):
    hits_df = hits_df.copy()
    hits_df["sequence"] = ""

    # Agrupa por genoma para leer cada .fna una sola vez
    for genome_dir, group in hits_df.groupby("genome_dir"):
        fna = glob.glob(os.path.join(genome_dir, "*.fna"))[0]
        contigs = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))

        for idx, row in group.iterrows():
            # Coordenadas BLAST: 1-based e inclusivas -> slicing de Python
            seq = contigs[row["seq_id"]].seq[row["start"] - 1 : row["end"]]
            if row["strand"] == "-":
                seq = seq.reverse_complement()  # orienta como la query
            hits_df.at[idx, "sequence"] = str(seq)

    return hits_df

def process_hits(tsv_path, dir_genomes, multifasta):
    hits_df = process_blast_tsv(tsv_path)  # con start, end, strand
    contig2genome = contig_to_genome(dir_genomes)

    # Solo las dos comprobaciones que evitan borrar todo por error
    if hits_df.empty:
        raise ValueError("BLAST no devolvió ningún hit: no borro nada.")
    if not hits_df["seq_id"].isin(contig2genome).all():
        raise ValueError("Hay seq_id sin genoma asignado: no borro nada.")

    hits_df["genome_dir"] = hits_df["seq_id"].map(contig2genome)
    hits_df["genome"] = hits_df["genome_dir"].map(os.path.basename)

    genomes_with_hit = set(hits_df["genome_dir"])
    to_delete = sorted(set(contig2genome.values()) - genomes_with_hit)

    for g in to_delete:
        shutil.rmtree(g)

    os.remove(multifasta)
    df = hits_df[["genome", "genome_dir", "seq_id", "start", "end", "strand"]]
    final_df = extract_hit_sequences(df)
    return final_df

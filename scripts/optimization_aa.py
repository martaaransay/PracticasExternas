import os
import subprocess
import time
import shlex
import glob
from pipeline_inicial.T3_2_integronfinder import Integron_Finder
from pipeline_inicial.T3_1_integronfiltering import integron_filtering
# FIXME: SOLUCIONAR DONDE ESTÁ LA FUNCIÓN DE LIST FASTA FILES 


def list_fasta_files(directory, num_directories = 1, extension = "*.fna"):
    """
    Gets the paths of all FASTA files (.fna) stored within the subfolders of a directory.

    Args
    -------
    directory (str): Path of the main directory containing the genome subfolders.

    Returns
    --------
    genome_paths (list): A list with the full paths to each found fasta (.fna) file.
    """
    # Check genome files
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"The directory {directory} was not found.")
    wildcard = ""
    for n in range(0, num_directories):
        if n>0: wildcard += "/"
        wildcard += "*"
        
    # Set the genome file path pattern
    genomes_path = os.path.join(directory, wildcard, extension)
    # Use glob to correctly intepret the wildcard (*)
    genomes_path_glob = sorted(glob.glob(genomes_path))
    
    return genomes_path_glob

if __name__ == "__main__":
    flag_integron_filtering = True
    
    
    n_genomes = 1500
    directory = f"data/{str(n_genomes)}genomes"
    list_files = list_fasta_files(directory, num_directories=5)
    start = time.time()
    dir_IF2 = f"results/{str(n_genomes)}genomes/IF2"
    dir_IntFilt = f"results/{str(n_genomes)}genomes/IntFilt"
    hits = {}
    os.makedirs(dir, exist_ok=True)
    for each_file in list_files:
        acc = os.path.splitext(os.path.basename(each_file))[0]
        if flag_integron_filtering:
            file = integron_filtering(genome_path = each_file,
                                   outdir=dir_IntFilt)
            if file:
                each_file = file
                acc = os.path.splitext(os.path.basename(each_file))[0]

        gbk_files = Integron_Finder(genome_path = each_file,
                        genome_acc = acc,
                        outdir = dir_IF2)
        hits[acc] = len(gbk_files) # cada elementos on los hits x genoma

    general_hits = sum(hits.values())
    unique_hits = 0
    for n in hits.values():
        if n > 0:
            unique_hits += 1
    
    total_time = time.time() - start #
    print(f"Total analyzed : {len(hits.keys(()))}")
    print(f"Total hits = {general_hits}")
    print(f"Total unique hits = {unique_hits}")

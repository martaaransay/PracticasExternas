# Import libraries 
import glob
import os
import shlex
import subprocess
import shutil
from pipeline_inicial.T2_downloadgenomes import datasets_flags, download_genomes

# --------Definition of functions
#------------------------ Function to download genomes for specific species
def download_species_genomes(species, outdir, target_num, flags):
    """
    Uses the download_genomes function to download genomes for each species in the species dictionary.
    
    Args
    ------
    species (dict): A dictionary where keys are short names and values are long names of species.
    outdir (str): The output directory where genomes will be downloaded.
    target_num (int): The target number of genomes to download for each species.
    flags (dict): Flags for the genome download.

    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok = True)
    # Iterates over the dictionary
    for short_name, long_name in species.items():
        # Generates the path to save the genomes
        path = os.path.join(outdir, short_name, f"{short_name}.zip")
        # Executes download function
        download_genomes(source = "taxon",
                        name = long_name, # Species to download
                        quality_filter_flag = True,
                        filename = path,
                        flags = flags,
                        target_num = target_num # number of genomes to download
                        )
        
    return

#------------------------ Function to generate a multifasta file
def generate_multifasta(path_to_genomes, outfile):
    """
    Generates a multifasta file with the concatenation of the fasta files of a specific path.

    Args
    ------
    path_to_genomes (str): The path to the genomes to be concatenated.
    outfile (str): The output file where the multifasta will be saved.
    """
    # Process the path to genomes to evaluate the wildcard correctly
    files_to_db = sorted(glob.glob(path_to_genomes))
    # Writes the multifasta file
    with open(outfile, "w") as out:
        for file in files_to_db: # For each file
            with open(file, "r") as each_file:
                out.write(each_file.read()) # Writes each fasta in the multifasta
                out.write("\n") # To avoid errors adds a new line between each fasta
    return


#------------------------ Function to generate the database from the multifasta file
def generate_db(multifasta, out_db):
    """
    Generates a BLAST database from a multifasta file.

    Args
    ------
    multifasta (str): Path to the multifasta file.
    out_db (str): Path where the database will be saved.
    """
    # The output directory is created to avoid errores
    os.makedirs("../data/genomes/db", exist_ok=True)
    # Command to execute makeblastdb
    makeblastdb_cmd = ["makeblastdb", "-in", multifasta, 
                       "-dbtype", "nucl", # nucleotides database
                       "-out", out_db, 
                       "-parse_seqids"] # To include the sequence identifiers in the database
    print(f"Executing: \n{shlex.join(makeblastdb_cmd)}")
    # Execute the command
    try:
        subprocess.run(makeblastdb_cmd, 
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing makeblastdb:\n{e.stderr}")
        return
    return

#------------------------ Function to remove all the files related to the database
def clean_db(outdir):
    """
    Removes all the files related to the database, including multifasta files, genomes, and the database itself.

    Args
    ------
    outdir (str): The output directory where the files to be removed are located.
    """
    shutil.rmtree(outdir) # Removes the entire directory and its contents

    return

# -------- Main code
if __name__ == "__main__":

    clean_flag = False
    if clean_flag:
            clean_db(outdir = "data/genomes")

    # Downloading genomes
    # Flags:
    flags = datasets_flags(assembly_level="complete,chromosome,contig",
                           exclude_atypical=True,
                           mag=False)
    n_genomes = 10
    # species to download
    names = {"Ecoli" : "Escherichia coli", 
             "Klebsiella" : "Klebsiella pneumoniae",
             "Pseudomonas" : "Pseudomonas aeruginosa",
             "Acinetobacter" : "Acinetobacter baumannii",
             "Enterobacter" : "Enterobacter"}

        
    download_species_genomes(species = names, 
                   outdir = f"data/{str(n_genomes)}genomes", 
                   target_num = n_genomes, 
                   flags = flags)

    
    generate_multifasta(path_to_genomes = f"data/{str(n_genomes)}genomes/*/*/ncbi_dataset/data/*/*.fna",
                        outfile = f"data/{str(n_genomes)}genomes/multifasta_genomes.fasta")
    generate_db(multifasta = f"data/{str(n_genomes)}genomes/multifasta_genomes.fasta",
                out_db = f"data/{str(n_genomes)}genomes/db/db")
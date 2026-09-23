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


# FIXME _ FALTA COMENTAR A PARTIR DE AQUÍ:
def generate_db(multifasta, out_db):
    # Generate db
    os.makedirs("../data/genomes/db", exist_ok=True)
    makeblastdb_cmd = ["makeblastdb", "-in", multifasta, "-dbtype", "nucl", 
                       "-out", out_db, "-parse_seqids"]
    print(f"Executing: \n{shlex.join(makeblastdb_cmd)}")
    try:
        subprocess.run(makeblastdb_cmd, 
                    check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing makeblastdb:\n{e.stderr}")
        return

def clean_db(outdir):
    #Removes multifasta + genomes + db
    shutil.rmtree(outdir)
    print("Everything removed")
    return

if __name__ == "__main__":
    clean_flag = False
    if clean_flag:
            clean_db(outdir = "data/genomes")

    # Downloading genomes
    flags = datasets_flags(assembly_level="complete,chromosome,contig",
                           exclude_atypical=True,
                           mag=False)
    
    names = {"Ecoli" : "Escherichia coli", 
             "Klebsiella" : "Klebsiella pneumoniae",
             "Pseudomonas" : "Pseudomonas aeruginosa",
             "Acinetobacter" : "Acinetobacter baumannii",
             "Enterobacter" : "Enterobacter spp."}
    
    download_species_genomes(species = names, 
                   outdir = "data/5000genomes", 
                   target_num = 5000, 
                   flags = flags)
    # Executed once
    generate_multifasta(path_to_genomes = "data/5000genomes/*/*/ncbi_dataset/data/*/*.fna",
                        outfile = "data/5000genomes/multifasta_genomes.fasta")
    generate_db(multifasta = "data/5000genomes/multifasta_genomes.fasta",
                out_db = "data/5000genomes/db/db")
    
    

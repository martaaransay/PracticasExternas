# Script with useful functions
import os, glob, re

def list_fasta_files(directory, depth = 1, extension = "*.fna"):
    """
    Gets the paths of all FASTA files (.fna) stored within the subfolders of a directory.

    Args
    -------
    directory (str): Path of the main directory containing the genome subfolders.
    depth (int, optional): Number of subfolder levels. 
                           Defaults to 1.
    extension (str, optional): Glob pattern used to match the fasta files. 
                               Defaults to "*.fna".

    Return
    --------
    genome_paths (list): A list with the full paths to each found fasta (.fna) file.
    """
    # Check genome files
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"The directory {directory} was not found.")
    wildcard = ""
    for n in range(0, depth):
        if n>0: wildcard += "/"
        wildcard += "*"
        
    # Set the genome file path pattern
    genomes_path = os.path.join(directory, wildcard, extension)
    # Use glob to correctly intepret the wildcard (*)
    genomes_path_glob = sorted(glob.glob(genomes_path))
    
    return genomes_path_glob


def contig_to_genome(dir_genomes):
    """
    Creates a mapping dictionary between contig IDs and their parent genome directories.
    
    Args
    --------
    dir_genomes (str): The root directory containing the downloaded NCBI datasets.

    Return
    --------
    contig2genome (dict): A dictionary where keys are contig IDs and values are the paths 
                          to the directory containing the genome FASTA file.
    """
    contig2genome = {} # Dictionary to store results
    # Pattern of each genome path
    pattern = os.path.join(dir_genomes, "*", "*", "ncbi_dataset", "data", "*", "*.fna")
    
    for fna in glob.glob(pattern): # Iterates over fasta fles
        genome_path = os.path.dirname(fna)  # Path
        with open(fna) as f:
            for line in f:
                if line.startswith(">"): # For the identifier,
                    id = line[1:] # Remove the '>'
                    id = id.split()[0] # Take the string before a blank space
                    # Map the id to the path
                    contig2genome[id] = genome_path
    return contig2genome

#------------------------ Function to normalize the accession between results of methods
def genome_key(name):
    """
    Extracts and normalizes the NCBI genome accession number from a given string.
        
    Args
    ------
    name (str): A genome path (without extension).

    Return
    ------
    accession (str): The extracted accession number if the pattern is found. 
                     If no match is found, it returns the original string.
    """
    pattern = r"GC[FA]_(\d+\.\d+)" # r"  to do not escape \d and \.
    match = re.search(pattern, str(name)) # Search the pattern into the path
    if match:
        accession = match.group(1) # Extracts the match
    else:
        accession = str(name)
    return accession

#------------------------ Function to deduplicate files
def unique_genome_files(files):
    """
    Keeps one file per assembly (GCF and GCA versions of a genome are the same).
    
    Args
    ------
    files (list): List of files to deduplicate.

    Return
    ------
    unique (dict): Dictionary with the genome accession as the key and the path
                   to the file as the value
    """
    unique = {}
    for file in sorted(files):
        unique.setdefault(genome_key(file), file) # Keeps the first one
    return unique 

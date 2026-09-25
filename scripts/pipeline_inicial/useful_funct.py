# Script with useful functions
import os, glob
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

    Returns
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

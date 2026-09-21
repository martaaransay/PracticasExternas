# Import libraries 
import subprocess # to execute cmd
import shlex # to check what is sent to the console
import json # to parse .json correctly
import os 

# --------Definition of functions
#------------------------ Function to define flags for genome download (datasets download)
#FIXME -> He dejado todas las opciones posibles que acepta datasets, pero pdoemos reducir las opciones¿¿¿¿ hasta cuanto vamos a considerar???
def datasets_flags(assembly_source = None, # RefSeq | GenBank | None =all
                   assembly_level = None, # comma-separated (string): chromosome | complete | contig | scaffold | None =all
                   annotated = None, # To exclude annotated genomes  (≠ None)
                   exclude_atypical = None, # Excludes atypical assemblies (≠ None) 
                   exclude_multi_isolate = None, # Excludes assemblies from multi-isolate projects  (≠ None)
                   mag = None, # If True or "only", only uses MAGS; if false or "exclude", excludes them
                   ):
    """
    Builds and validates a list of flags for the NCBI datasets download command.

    Args
    --------
    assembly_source (str): Assembly source. 
                           Allowed values: 'RefSeq' or 'GenBank' (eg. assembly_source='RefSeq').
    assembly_level (str): Assembly level - comma-separated without blank.
                          Allowed values: 'chromosome', 'complete', 'contig', 'scaffold' (eg. assembly_source='chromosome,complete').
    annotated (bool): If evaluated to ≠ False/None, includes annotated genomes ('--annotated').
    exclude_atypical (bool): If evaluated to ≠ False/None, excludes atypical assemblies ('--exclude-atypical').
    exclude_multi_isolate (bool): If evaluated to ≠ False/None, excludes multi-isolate projects ('--exclude-multi-isolate').
    mag (bool | str): Filter for metagenomic genomes. 
                      Allowed values: True/'only', False/'exclude' or 'all' (eg. mag='only').

    Return
    --------
    flags (list): List of strings with the arguments to be passed to subprocess.
    """
    
    flags = [] # List of flags to fill

    if assembly_source: # If an argument was provided, adds it to the flags list
        # datasets has a checkpoint control in case the argument is not valid
        flags += ["--assembly-source", assembly_source] # Add the argument to the list
        
    if assembly_level: 
        # datasets has a checkpoint control in case the argument is not valid
        flags += ["--assembly-level", assembly_level]

    if annotated:
        flags.append("--annotated")

    if exclude_atypical:
        flags.append("--exclude-atypical")

    if exclude_multi_isolate:
        flags.append("--exclude-multi-isolate")

    if mag is not None: 
        if mag is True:
            flags += ["--mag", "only"]
        elif mag is False:
            flags += ["--mag", "exclude"]
        else:
            flags += ["--mag", mag]

    return flags


#------------------------ Function to define taxon/accesion 
def taxon_acc_flag(source, # taxon | accession | inputfile
                   name # List of values (source="taxon" | source="accession")\
                        # or file name (source="inputfile")
                   ):
    """
    Builds the argument list (taxon or accession info) for the NCBI datasets download command.

    Args
    --------
    source (str): Input type for the search. 
                  Allowed values: 'taxon', 'accession' or 'inputfile'.
    name (str): Values to search (comma-separated string) or file path if the source is 'inputfile'.

    Return
    --------
    source_flags (list): List of strings with the unpacked source arguments for datasets.
    """
    
    if source not in ["taxon", "accession", "inputfile"]: # Control of valid arguments (inputfile to internal control is added)
        raise ValueError(f"One of the valid options was not provided: taxon | accession | inputfile \
                        \n (provided: {source})")
    
    if source == "inputfile": # File
        if not os.path.isfile(name): # Check in case the file does not exist
            raise FileNotFoundError(f"The file {name} was not found.")
            return
        flag = ["accession", "--inputfile", name] # Flags to add with the name of the file
    
    else:
        names_flag = [] 
        for each_name in name.split(","): # Single element or several comma-separated
            names_flag.append(each_name)

        flag = [source, *names_flag] # * unpacks, so it does not generate a nested list
    return flag

#------------------------ Function to filter genomes
def quality_filter_genomes(source, name, # Arguments for the taxon_acc_flag(source, name) function
                           filename, # Filename to keep the same directory
                           min_n50 = 2000, 
                           max_contamination = 5, 
                           min_completeness = 90): 
    """
    Queries genome metadata in NCBI datasets, filters by quality metrics, and saves valid accessions in a text file.

    Args
    ------------
    source (str): Source type for the query.
    name (str): Value or identifier associated with the source.
    filename (str): Name of the resulting .zip file. 
    n50_threshold (int, optional): Minimum acceptable value for the scaffold N50. 
                                   Default is 2000.
    max_contamination (float | int, optional): Maximum tolerated contamination percentage according to CheckM. 
                                               Default is 5%.
    min_completeness (float | int, optional): Minimum required completeness percentage according to CheckM. 
                                                    Default is 90%.

    Return
    --------------
    accessions_tuple (tuple): Tuple with the generated file name if there are genomes that pass the filters
                              None if there are no valid results.
    """
    # To get the JSON with the genome data 
    summary_cmd = ["datasets", "summary", "genome", 
                    *taxon_acc_flag(source, name), "--as-json-lines"]
    
    print(f"Executing: \n{shlex.join(summary_cmd)}")
    #FIXME: tomar decisión de como gestionar las excepciones!
    try:
        summary_json = subprocess.run(summary_cmd, 
                                      capture_output = True, # Captures stdout and stderr
                                      text = True, # stdout and stderr files in text mode
                                      check = True) # To capture if an error occurs (CalledProcessError exception)
    except subprocess.CalledProcessError as e: # Captures the exception
        print(f"Error executing datasets summary:\n{e.stderr}")
        return

    valid_accs = [] # List to store valid accessions 
    summary_text = summary_json.stdout.splitlines() # .json obtained (separated by lines)

    for line in summary_text:
        if not line.strip(): continue # Skips empty lines
        datos = json.loads(line) # Deserialization
        acc = datos.get("accession") # Get the accession
        
        # Quality metrics:
        #------------------- N50
        assembly_stats = datos.get("assembly_stats", {}) # Safely get stats (if not present, empty dict)
        n50 = assembly_stats.get("scaffold_n50", 0)
        
        #------------------- CheckM info
        checkm_info = datos.get("checkm_info", {})
        # If there is no checkm info, the default values are set
        # Contamination 
        contamination = checkm_info.get("contamination", 5) 
        # Completeness
        completeness = checkm_info.get("completeness", 90) 

        # Check if it meets the selected filters 
        if n50 >= min_n50 and contamination <= max_contamination and completeness >= min_completeness:
            valid_accs.append(acc)

    if not valid_accs:
        print("No genome meets the quality criteria.")
        return
    
    accs_file_name = f"accessions_{os.path.basename(filename).split('.')[0]}.txt" # Name of the file
    accs_file_path = os.path.join(os.path.dirname(filename), accs_file_name)
    print(accs_file_path)
    
    with open(accs_file_path, "w") as f: # Create a file with the valid accs (one per line)
        for each_acc in valid_accs:
            f.write(f"{each_acc}\n")
        print(f"The accessions to be analyzed have been saved in the file '{accs_file_path}'.")

    accessions_tuple = ("inputfile", accs_file_path)

    return accessions_tuple
##------------------------ Function to unzip and rehydrate the genomes
def unzip_rehydrate(filename = "ncbi_dataset.zip"): 
    """
    Unzips an NCBI ZIP file and executes the datasets rehydrate command.

    Args
    -------
    filename (str, optional): Path to the ZIP file to decompress. 
                              Default is "ncbi_dataset.zip".
    """
    if not os.path.isfile(filename): # Check that the file exists
        raise FileNotFoundError(f"The file {filename} was not found")

    ############# ------------------- unzip
    # If there' is no path in the filename, returns "."
    out_dir = os.path.dirname(filename) or "." 
    new_file = os.path.basename(filename)
    new_folder = os.path.splitext(new_file)[0]
    target_dir = os.path.join(out_dir, new_folder)

    unzip_cmd = ["unzip",
                 "-q", # Quiet
                 "-o", # Overwrite files to avoid errors
                 filename,
                 "-d", target_dir] # Directory where it is unzipped
    
    print(f"Executing: \n{shlex.join(unzip_cmd)}")
    try:
        subprocess.run(unzip_cmd, 
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing unzip:\n{e.stderr}")
        return

    ############ ------------------- datasets rehydrate
    rehydrate_cmd = ["datasets", "rehydrate", 
                     "--directory", target_dir]
    print(f"Executing: \n{shlex.join(rehydrate_cmd)}")
    try:
        subprocess.run(rehydrate_cmd, 
                        check=True, capture_output=True, text=True) 
    except subprocess.CalledProcessError as e:
        print(f"Error executing datasets rehydrate:\n{e.stderr}")
        return

    return


##------------------------ Function to download genomes
def download_genomes(source, name, # Arguments for the taxon_acc_flag(source, name) function
                     quality_filter_flag = False, # Flag to filter by quality metrics is needed
                     filename = "ncbi_dataset.zip", # Name of the file that will contain the genomes 
                                                  # Can have a path to specify a target directory
                     flags = None): # Flags for the datasets command
    """
    Executes the NCBI datasets command to download genomes (--dehydrated), allowing prior quality filtering.

    Args
    --------
    source (str): Source type for the query.
    name (str): Value or identifier associated with the source.
    quality_filter_flag (bool, optional): If evaluates to True, executes quality_filter_genomes() function.
                                          Default is False.
    filename (str, optional): Name of the resulting .zip file.  
                              Default is "ncbi_dataset.zip".
    flags (list, optional): List of additional flags to add to the 'datasets download' command. 
                            Default is None.
    """
    # Check that the directories exist to prevent errors
    out_dir = os.path.dirname(filename)
    if out_dir: # Prevents errors if the path has no folders 
        os.makedirs(out_dir, exist_ok=True) 

    if quality_filter_flag:
        filter_result = quality_filter_genomes(source, name, filename) # Executes the filtering function
        if not filter_result: # If there is no accession, the datasets command cannot be used
            # The created directories are removed:
            if out_dir and os.path.exists(out_dir):
                if not os.listdir(out_dir):
                    # Removes empty directories (parent directories too) 
                    os.removedirs(out_dir)
            return
        source, name = filter_result # Unpack tuple

    # If no flags were added, empty list (control)
    if flags is None:
        flags = []

    ########### -------------------- datasets download
    download_cmd = ["datasets", "download", "genome", *taxon_acc_flag(source, name), 
                    "--dehydrated", "--filename", filename, *flags]
    print(f"Executing: \n{shlex.join(download_cmd)}")

    try:
        subprocess.run(download_cmd, 
                       check=True, capture_output=True, text=True) 
    # Capture the exception
    except subprocess.CalledProcessError as e:
        print(f"Error executing datasets download:\n{e.stderr}")
        # The created directories are removed (provided they are completely empty)
        if out_dir and os.path.exists(out_dir):
            if not os.listdir(out_dir): 
                os.removedirs(out_dir) # Removes empty directories (parent directories too)
        return None
    unzip_rehydrate(filename)


# -------- Main code

if __name__ == "__main__":
    # Code for Klebsiella pneumoniae without quality filtering
    flags = datasets_flags(assembly_source = "RefSeq", 
                        assembly_level = "complete,chromosome",
                        exclude_atypical = True, 
                        mag = False)

    file = "results/T2_check/kpn_without_qfiltering.zip"

    # check = download_genomes(source = "taxon", 
    #                          name = "Klebsiella pneumoniae", 
    #                          filename = file, 
    #                          flags = flags)

    # if check:
    #     unzip_rehydrate(filename = file)

    # Code for Klebsiella pneumoniae WITH quality filtering
    filtered_file = "results/T2_check/kpn_filtered.zip"

    # filtered_check = download_genomes(source = "taxon", 
    #                                   name = "Klebsiella pneumoniae", 
    #                                   quality_filter_flag = True,
    #                                   filename = filtered_file, 
    #                                   flags = flags)

    # if filtered_check:
    #     unzip_rehydrate(filename = filtered_file)


# FIXME -> no salen las barras de progreso al ejecutar este script, pero si quitamos el argumento (capture_output=True)
#          en subprocess.run, sale -> solo que ya no se captura el error así (habría que eliminar esa parte de:
#          (except subprocess.CalledProcessError as e:)
#          PERO OJO -> SOLO EN LAS QUE HAGA FALTA (EJ. en unzip no!)

# ###########PRUEBAS 16/09  -(habrá que eliminarlo):
# # Para probar si IntegronFinder2 solo identifica P2
# file01 = "results/pruebas_IntegronFinderCheck/EC.zip"
# flags01 = definir_flags(assembly_source="RefSeq", 
#                       assembly_level=["complete", "chromosome"],
#                       exclude_atypical=True, 
#                       mag=False)

# check01 = descargar_genomas(origen="taxon", 
#                           nombre="Escherichia coli", 
#                           filename=file01, 
#                           flags=flags01)

####PRUEBAS CLUSTERIZADO 17/09

if __name__ == "__main__":
    # Code for Klebsiella pneumoniae without quality filtering
    flags_general = datasets_flags()

    file_general = "results/T2_check/Klebsiella.zip"

    download_genomes(source = "taxon", 
                     name = "Klebsiella pneumoniae", 
                     filename = file_general, 
                     flags = flags_general)

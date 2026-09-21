# Import libraries 
import os
import re
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# --------Definition of functions
#------------------------ Function to extract information about integrons from a gbk from IntegronFinder2
def extract_integron_info_gbk(genome_acc, gbk_files):
    """
    Extracts integron information from .gbk files.

    Args:
    ------
        genome_acc (str): The accession number of the genome.
        gbk_files (list): A list containing the file paths to the .gbk files.

    Return:
    ------
        integron_info (list of dict): A list of dictionaries containing the extracted information. 
            For CALIN integrons, the dictionary contains:
                - "acc" (str): Genome accession number.
                - "calin" (bool): True.
                - "integron_id" (str): The ID of the integron.
                - "locus_id" (str): The replicon/locus ID.
            For non-CALIN integrons, it also includes:
                - "calin" (bool): False.
                - "start" (int): Start position of the integrase feature.
                - "end" (int): End position of the integrase feature.
                - "sequence" (Bio.Seq.Seq): The extracted sequence of the integrase.
    """
    # List to store integron information
    integron_info = []
    # Iterate through the .gbk files
    for gbk_file in gbk_files: # (more than one file can be generated from 1 genome)
        record = SeqIO.read(gbk_file, "genbank") # Read .gbk file
        locus_id = record.id # Store replicon id

        flag_calin = False  # Flag to control CALIN
   
        for feature in record.features: # Iterate through features

            if feature.type == "integron":
                # Get integron id (in case no value, empty list)
                integron_id = feature.qualifiers.get("integron_id", [""])[0]
                # Get integron type
                type_integron = feature.qualifiers.get("integron_type", [""])[0]

                if type_integron == "CALIN":
                    flag_calin = True # Changes flag
                    # Store info about calin-integron
                    integron_info.append({"acc": genome_acc,
                                          "calin": True,
                                          "integron_id": integron_id,
                                          "locus_id": locus_id
                    })

            if not flag_calin:  # no CALIN (integron_type -> complete/incomplete)      
                if feature.type == "integrase": # Pc is in the integrase sequence

                    start = int(feature.location.start) # Start of Int
                    end = int(feature.location.end) # End of Int

                    seq = record.seq[start:end] # Sequence of Int

                    # Store info about integron
                    integron_info.append({"acc": genome_acc,
                                          "calin": False,
                                          "integron_id": integron_id,
                                          "locus_id": locus_id,
                                          "start": start,
                                          "end": end,
                                          "sequence": seq
                        }) 
    
    return integron_info

# Regular Expressions of the pc variants
pc_variants_regex = {
    "PcS": ("TTGACA[ACGT]{17}TAAACT"),
    "PcW": ("TGGACA[ACGT]{17}TAAGCT"),
    "PcWTGN-10": ("TGGACA[ACGT]{14}TG[ACGT]TAAGCT"),
    "PcH1": ("TGGACA[ACGT]{17}TAAACT"),
    "PcH1TGN-10": ("TGGACA[ACGT]{14}TG[ACGT]TAAACT"),
    "PcH2": ("TTGACA[ACGT]{17}TAAGCT"),
    "PcH2TGN-10": ("TTGACA[ACGT]{14}TG[ACGT]TAAGCT"),
    "TGGGCA-N14-TGn-TAAGCT": ("TGGGCA[ACGT]{14}TG[ACGT]TAAGCT"),
    "PcSS": ("TTGATA[ACGT]{17}TAAACT"),
    "Pcln42": ("TTGGCA[ACGT]{17}TAAACT"),
    "Pcln116": ("TTGACA[ACGT]{17}TGAACT"),
    "PcPUO": ("TCGACA[ACGT]{17}TAAACT"),
    "P2": ("TTGTTA[ACGT]{17}TACAGT"),
    "P2m1":("TTGTTA[ACGT]{17}GACAGT"),
    "P2m2": ("TTGTTA[ACGT]{17}TACACA")
}

#------------------------ Function to search and classify Pc in Integrase sequences
def classify_pc_regex(sequence, variants = pc_variants_regex):
    """
    Searches and classifies Pc promoter variants within an integrase sequence using regular expressions. 
    
    Args:
    ------
    sequence (str): The nucleotide sequence of the integrase to be analyzed.
    variants (dict, optional): A dictionary containing:
                                -key: names of the Pc variant.
                                -value: their corresponding regular expression patterns. 
                               Defaults to pc_variants_regex.

    Return:
    ------
    results (dict): A dictionary containing the classification results of the matched Pc variant. 
                    Returns None if no match is found. 
    """
    # As the Pc is in the opposite direction of Integrase, the reverse complement is analyzed
    seq = str(Seq(sequence).reverse_complement())
    # Initialize a variable to store results
    results = None 
    # Get the name of the Pc and the RE of the sequence (.items() returns key and value)
    for name_pc_variant, pc_sequence in variants.items():
        # Find the pc sequence into the RC seq of the Integrase
        matches = re.finditer(pc_sequence, seq) 
        for match in matches: # For each match, stores the name, the positions and the sequence
            results = {
                "Pc_variant": name_pc_variant,
                "start": match.start(),
                "end": match.end(),
                "sequence": match.group()
                }
    return results

#------------------------ Function to store information from dictionaries into a csv file
def info_to_csv(dict_info, path_out_file):
    """
    Appends the information from a dictionary into a CSV file. 
    
    Args:
    ------
    dict_info (dict): A dictionary containing the data to be stored. 
    path_out_file (str): The file path where the CSV file will be saved.
    """
    # Open the file in appending mode
    with open(path_out_file, "a", newline="") as out_csv:
        csvwriter = csv.writer(out_csv) # Writer
        # If the file is new, writes as header the keys of the dictionary
        if os.path.getsize(path_out_file) == 0:
            csvwriter.writerow(dict_info.keys())
        # Writes the values of the dictionary
        csvwriter.writerow(dict_info.values())
    return

#------------------------ Function to merge general and Pc information into a single dictionary
def adjust_Pc_info(general_info, # From extract_integron_info_gbk() function
                   Pc_info, # From classify_pc_regex() function
                   path_out_file, flag_Pc_seq = None): # Flag to store the sequence of the Pc 
    """
    Merges general integron information and Pc variant classification details into a single 
    dictionary and appends the resulting data to a CSV file.

    Args:
    ------
    general_info (dict): A dictionary containing general integron information.
                         It must include the keys "acc", "integron_id", and "locus_id".
    Pc_info (dict): A dictionary containing the classification results of the Pc variant. 
                    It must include the keys "start", "end", and "Pc_variant", and 
                    "sequence" if the sequence is to be stored.
    path_out_file (str): The file path where the merged information will be saved as a CSV.
    flag_Pc_seq (bool, optional): A flag indicating whether to include the exact Pc sequence. 
                             Defaults to None.
    """
    # Extract values from general_info dictionary
    acc_genome = general_info["acc"]
    integron_id = general_info["integron_id"]
    locus_id = general_info["locus_id"]

    # Extract values from Pc_info dictionary
    start = Pc_info["start"]
    end = Pc_info["end"]
    Pc_variant = Pc_info["Pc_variant"]

    # Final dictionary with all the information
    final_results = {"acc_genome" : acc_genome,
                     "integron_id" : integron_id,
                     "locus_id" : locus_id,
                     "start_Pc": start,
                     "end_Pc": end,
                     "Pc_variant": Pc_variant}
    # If flag is True, the sequence of the Pc is added
    if flag_Pc_seq:
        final_results["sequence"] = Pc_info["sequence"]

    # Stores the final dictionary into the csv file
    info_to_csv(final_results, path_out_file)
    return 

#------------------------ Function to sort general information into a dictionary
def adjust_calin_info(general_info, # From extract_integron_info_gbk() function
                      path_out_file):
    """
    Formats and extracts specific general information for CALIN integrons 
    into a new dictionary, and appends the resulting data to a CSV file.

    Args:
    ------
    general_info (dict): A dictionary containing general integron information. 
                         It must include the keys "acc_genome", "integron_id", and "locus_id".
    path_out_file (str): The file path where the extracted information will be saved as a CSV.
    """
    
    # Extract values from general_info dictionary
    acc_genome = general_info["acc"]
    integron_id = general_info["integron_id"]
    locus_id = general_info["locus_id"]

    # Final dictionary with all the information
    final_results = {"acc_genome" : acc_genome,
                     "integron_id" : integron_id,
                     "locus_id" : locus_id}

    # Stores the final dictionary into the csv file
    info_to_csv(final_results, path_out_file)

    return 


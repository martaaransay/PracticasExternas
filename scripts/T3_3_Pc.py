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
                - "sequence" (Bio.Seq.Seq): The extracted sequence of the integrase and its upstream region.
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


variantes_pc = {
    "PcS": ("TTGACA", 17, "TAAACT"),
    "PcW": ("TGGACA", 17, "TAAGCT"),
    "PcWTGN-10": ("TGGACA", 14, "TG", 1, "TAAGCT"),
    "PcH1": ("TGGACA", 17, "TAAACT"),
    "PcH1TGN-10": ("TGGACA", 14, "TG", 1, "TAAACT"),
    "PcH2": ("TTGACA", 17, "TAAGCT"),
    "PcH2TGN-10": ("TTGACA", 14, "TG", 1, "TAAGCT"),
    "TGGGCA-N14-TGn-TAAGCT": ("TGGGCA", 14, "TG", 1, "TAAGCT"),
    "PcSS": ("TTGATA", 17, "TAAACT"),
    "Pcln42": ("TTGGCA", 17, "TAAACT"),
    "Pcln116": ("TTGACA", 17, "TGAACT"),
    "PcPUO": ("TCGACA", 17, "TAAACT"),
    "P2": ("TTGTTA", 17, "TACAGT"),
    "P2m1":("TTGTTA", 17, "GACAGT"),
    "P2m2": ("TTGTTA", 17, "TACACA")
}


variantes_pc_regex = {
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


def clasificar_Pc_regex(seq, variantes = variantes_pc_regex):
    """
    """
    sequence = str(Seq(seq).reverse_complement())
    resultados = None
    for name_pc_variant, pc_sequence in variantes.items():
        matches = re.finditer(pc_sequence, sequence)
        for match in matches:
            resultados = {
                "variante": name_pc_variant,
                "start": match.start(),
                "end": match.end(),
                "secuencia": match.group()
            }

    return resultados

def info_to_csv(dict_info, path_out_file):
    
    with open(path_out_file, "a", newline="") as out_csv:
        csvwriter = csv.writer(out_csv)
        if os.path.getsize(path_out_file) == 0:
            #quiero que escriba las keys
            csvwriter.writerow(dict_info.keys())
        #quiero que escriba las keys
        csvwriter.writerow(dict_info.values())
    return

def adjust_Pc_info(general_info, Pc_info, path_out_file, Pc_seq = None):
    """docstring"""
    # from general info: (extract_pc_gbk...)

    acc_genoma = general_info["acc"]
    integron_id = general_info["integron_id"]
    locus_id = general_info["locus_id"]

    # from Pc_info (clasificar Pc regex)

    start = Pc_info["start"]
    end = Pc_info["end"]
    pc_variant = Pc_info["variante"]

    final_results = {"acc_genoma" : acc_genoma,
                     "integron_id" : integron_id,
                     "locus_id" : locus_id,
                     "start_Pc": start,
                     "end_Pc": end,
                     "Pc_variant": pc_variant}

    if Pc_seq:
        final_results["sequence"] = seq
    info_to_csv(final_results, path_out_file)
    return 

def adjust_calin_info(general_info, path_out_file):
    """docstring"""
    # from general info: (extract_pc_gbk...)

    acc_genoma = general_info["acc_genoma"]
    integron_id = general_info["integron_id"]
    locus_id = general_info["locus_id"]

    final_results = {"acc_genoma" : acc_genoma,
                     "integron_id" : integron_id,
                     "locus_id" : locus_id}
    info_to_csv(final_results, path_out_file)
    return 


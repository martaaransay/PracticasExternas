
import subprocess
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_Pc_region_gbk(genome_acc, gbk_files, upstream = 0):
    """docstring"""
    # Lista para almacenar la información de los promotores
    promotores_Pc = []

    for gbk_file in gbk_files: # Recorre los archivos .gbk
        # Lee el archivo .gbk
        record = SeqIO.read(gbk_file, "genbank")
        locus_id = record.id # Almacena el nombre del replicón
        # Control por si aparece antes Pc que integron y evitar error
        n_integron = 0
        for feature in record.features: # Recorre las features
            if feature.type == "integron":
                n_integron += 1
                flag_calin = False
                type_integron = feature.qualifiers.get("integron_type", [""])[0]
                if type_integron == "CALIN":
                    flag_calin = True
                    promotores_Pc.append({
                        "acc_genoma": genome_acc,
                        "calin": True,
                        "n_integron": n_integron,
                        "locus": locus_id
                    })
            if not flag_calin:              
                if feature.type == "integrase": # Almacena el nombre del integrón
                    # Si no existe, devuelve vacío para evitar errores; y coge el primer elemento de la lista
                    # integron_id = feature.qualifiers.get("integron_id",
                    #                                      [""])[0]
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    start_upstream = start - upstream
                    if start_upstream <= 0: start_upstream = 1
                    seq = record.seq[start_upstream:end]

                    promotores_Pc.append({
                            "acc_genoma": genome_acc,
                            "calin": False,
                            "n_integron": n_integron,
                            "locus": locus_id,
                            "start": start_upstream,
                            "end": end,
                            "secuencia": seq
                        }) 

    return promotores_Pc

# esta función solo se ejecuta una vez para generar la base de datos y ya
def generar_db(fasta_file = "data/Pc/Pc.fasta",
               db_out = "data/Pc/db/Pc_db"):
    bbdd_cmd = ["makeblastdb", "-in", fasta_file, "-dbtype", "nucl", "-out", db_out]

    try:
        subprocess.run(bbdd_cmd, capture_output=True, text=True, check=True)
    
    except subprocess.CalledProcessError as e:
        print(f"Error executing makeblastdb:\n{e.stderr}")
        return None
    return


def seq_to_fasta(seq, acc, out, line_length=80):

    """docstring"""
    
    record = SeqRecord(Seq(seq), id=acc, description="" )
    acc = acc + ".fasta"
    out_file = os.path.join(out, acc)
    SeqIO.write(record, out_file, "fasta")
    return out_file

def clasificar_Pc(seq_fasta, db="data/Pc/db/Pc_db"):
    dir_out = os.path.dirname(seq_fasta)
    name = os.path.splitext(os.path.basename(seq_fasta))[0]
    name_out = name + ".out"
    file_out = os.path.join(dir_out, name_out)
    blast_cmd = ["blastn", "-db", db, "-query", seq_fasta, "-out", file_out, "-outfmt", "7"]
    try:
        subprocess.run(blast_cmd, capture_output=True, text=True, check=True)
    
    except subprocess.CalledProcessError as e:
        print(f"Error executing blast:\n{e.stderr}")
        return None
    print("correcto")
    return


variantes_pc = {
    "PcS":  ("TTGACA", 17, "TAAACT"),
    "PcW":  ("TGGACA", 17, "TAAGCT"),
    "PcH1": ("TGGACA", 17, "TAAACT"),
    "PcH2": ("TTGACA", 17, "TAAGCT"),
    "P2":   ("TTGTTA", 17, "TACAGT"),
}


def clasificar_Pc_regex(seq):
    """
    """
    seq = str(seq).upper()

    rc = str(Seq(seq).reverse_complement())


    resultados = []
    for strand_name, s in [("+", seq), ("-", rc)]:
        for nombre, (box35, spacer_len, box10) in variantes_pc.items():
            patron = f"{box35}[ACGT]{{{spacer_len}}}{box10}"
            for match in re.finditer(patron, s):
                resultados.append({
                    "variante": nombre,
                    "strand": strand_name,
                    "start": match.start(),
                    "end": match.end(),
                    "secuencia": match.group()
                })
    return resultados


import subprocess
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

def construir_patron(partes):
    box35 = partes[0]
    box10 = partes[-1]
    medio = partes[1:-1]

    patron = box35
    for parte in medio:
        if isinstance(parte, int):
            patron += f"[ACGT]{parte}"
            print(patron)
        else:
            patron += parte
    patron += box10

    return patron


def clasificar_Pc_regex(seq):
    """
    """
    seq = str(seq).upper()
    rc = str(Seq(seq).reverse_complement())

    resultados = []
    for strand_name, s in [("+", seq), ("-", rc)]:
        for nombre, partes in variantes_pc.items():
            patron = construir_patron(partes)
            for match in re.finditer(patron, s):
                resultados.append({
                    "variante": nombre,
                    "strand": strand_name,
                    "start": match.start(),
                    "end": match.end(),
                    "secuencia": match.group()
                })
    return resultados


clasificar_Pc_regex(Seq("CASDCFSADFSDFASFD"))
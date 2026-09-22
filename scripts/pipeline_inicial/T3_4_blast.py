#FIXME: CODIGO PARA EJECUTAR LA BUSQUEDA DE BLAST, PERO TRAS CAMBIAR VARIOS PARÁMETROS NO FUNCIONA BIEN POR LA CORTA LOGNITUD DE LOS PC
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

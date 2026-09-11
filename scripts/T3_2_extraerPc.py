import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def extraer_secuencia_Pc(ruta_genoma, acc_genoma, archivo_integrons):
    df = pd.read_csv(archivo_integrons, 
                     sep="\t",
                     comment="#")
    #    FIXME -> SE PUEDEN ELIMINAR SIN PORBLEMA OTRAS COLUMNAS Q NO SEAN RELEVANTES...
    promotores= df[(df["type_elt"] == "Promoter") & \
                   (df["element"].str.startswith("Pc"))].copy()
    promotores.insert(0, "ID_genome", acc_genoma)

    genoma_dict = SeqIO.index(ruta_genoma, "fasta")
    
    secuencias = []
    for index, row in promotores.iterrows():
        id_replicon = row['ID_replicon']
        
        # Verificar que el cromosoma/plásmido esté en el archivo FASTA
        if id_replicon not in genoma_dict:
            print(f"Advertencia: {id_replicon} no encontrado en el archivo FASTA. Saltando...")
            continue
            
        # Extraer variables de interés (asegurando formato entero)
        pos_beg = int(row["pos_beg"])
        pos_end = int(row["pos_end"])
        strand = int(row["strand"])
        
        # 3. Extraer la secuencia
        # Python recorta con [inicio:fin]. Restamos 1 a pos_beg para el índice 0.
        secuencia_cruda = genoma_dict[id_replicon].seq[pos_beg - 1 : pos_end]
        
        # 4. Ajustar según la dirección de la hebra (Strand)
        if strand == -1:
            # Obtiene el reverso complementario automáticamente
            secuencia_final = secuencia_cruda.reverse_complement()
        else:
            secuencia_final = secuencia_cruda
            
        # Guardar el resultado
        secuencias.append(str(secuencia_final))
        
    # 5. Añadir las secuencias como una nueva columna al DataFrame original
    df_resultado = promotores.copy()
    df_resultado["sequence"] = secuencias
    
    return df_resultado


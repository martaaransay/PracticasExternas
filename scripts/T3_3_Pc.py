from Bio import SeqIO

def extract_Pc_region_gbk(genome_acc, gbk_files, upstream = 100):
    """
    Extrae la información de los Pc a partir de archivos .gbk.

    Args
    ------
    genome_acc (str): Nombre de accesión del genoma.
    gbk_files (list): Lista de rutas de los archivos .gbk.

    Returns
    ------
    promotores_Pc (list[dict]): Lista de diccionarios.
                  Cada diccionario se corresponde a un Pc con las siguientes claves:
                    - "acc_genoma": (str) Accesión del genoma.
                    - "locus": (str) Identificador del replicón.
                    - "integron_id": (str) Identificador del integrón.
                    - "nombre": (str) Nombre del promotor.
                    - "start": (int) Coordenada de inicio.
                    - "end": (int) Coordenada de fin.
                    - "hebra": (int) Hebra.
                    - "secuencia": (Seq) Secuencia del promotor.
    """
    # Lista para almacenar la información de los promotores
    promotores_Pc = []

    for gbk_file in gbk_files: # Recorre los archivos .gbk
        # Lee el archivo .gbk
        record = SeqIO.read(gbk_file, "genbank")
        locus_id = record.id # Almacena el nombre del replicón
        # Control por si aparece antes Pc que integron y evitar error
        integron_id = ""

        for feature in record.features: # Recorre las features
            if feature.type == "integron": # Almacena el nombre del integrón
                # Si no existe, devuelve vacío para evitar errores; y coge el primer elemento de la lista
                integron_id = feature.qualifiers.get("integron_id",
                                                     [""])[0]
                start = int(feature.location.start)
                end = int(feature.location.end)
                start_upstream = start - upstream
                if start_upstream <= 0: start_upstream = 1
                seq = record.seq[start_upstream:end]
                promotores_Pc.append({
                        "acc_genoma": genome_acc,
                        "locus": locus_id,
                        "integron_id": integron_id,
                        "start": start,
                        "start_upstream": start_upstream,
                        "end": end,
                        "secuencia": seq
                    }) 

    return promotores_Pc

def clasificar_Pc():
    #FIXME: Aquí es donde hay que hacer la base de datos con BLAST
    return
from Bio import SeqIO

def extraer_Pc_gbk(acc_genoma, archivos_gbk_integronfinder):
    """
    Extrae la información de los Pc a partir de archivos .gbk.

    Args
    ------
    acc_genoma (str): Nombre de accesión del genoma.
    archivos_gbk_integronfinder (list): Lista de rutas de los archivos .gbk.

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

    for archivo_gbk in archivos_gbk_integronfinder: # Recorre los archivos .gbk
        # Lee el archivo .gbk
        record = SeqIO.read(archivo_gbk, "genbank")
        locus_id = record.id # Almacena el nombre del replicón
        # Control por si aparece antes Pc que integron y evitar error
        integron_id = ""

        for feature in record.features: # Recorre las features
            if feature.type == "integron": # Almacena el nombre del integrón
                # Si no existe, devuelve vacío para evitar errores; y coge el primer elemento de la lista
                integron_id = feature.qualifiers.get("integron_id",
                                                     [""])[0] 
            if feature.type == "Promoter": # Si es un promotor
                nombre = feature.qualifiers.get("Promoter", [""])[0] # Almacena su nombre
                if "Pc" in nombre: # Si es un Pc (NO promotor de la integrasa)
                    start = int(feature.location.start) # Coord inicio del promotor
                    end = int(feature.location.end) # Coord final del promotor
                    hebra = feature.location.strand # Hebra en la que se encuentra el promotor
                    seq = feature.location.extract(record.seq) # Secuencia del promotor
                    # Almacena la información extraída del promotor en la lista en formato de diccionario
                    promotores_Pc.append({
                        "acc_genoma": acc_genoma,
                        "locus": locus_id,
                        "integron_id": integron_id,
                        "nombre": nombre,
                        "start": start,
                        "end": end,
                        "hebra": hebra,
                        "secuencia": seq
                    })
    return promotores_Pc

def clasificar_Pc():
    return
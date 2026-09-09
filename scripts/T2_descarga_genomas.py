import subprocess
# subprocess permite ejecutar comandos en la terminal (ej. datasets)

###FUNCIÓN PARA DEFINIR EL FILTRADO DEL GENOMA

def definir_flags(
          assembly_source = None, # RefSeq | GenBank | None =todas
          assembly_level = None, # comma-separated: chromosome | complete | contig | scaffold | None =todas
          annotated = None, # Si no es None, excluye genomas anotados
          exclude_atypical = None, # Si no es None, excluye assemblies atípicos
          exclude_multi_isolate = None, # Si no es None, excluye assemblies de proyectos multi-isolate
          mag = None, # Si es True, solo usa mags; si es false, los excluye
          ):
    
    # Lista de las flags

    flags = []

    if assembly_source:
        flags += ["--assembly-source", assembly_source]
    
    if assembly_level:
        # Comma-separated:
        flags += ["--assembly-level", ",".join(assembly_level)]

    if annotated:
        flags.append("--annotated")

    if exclude_atypical:
        flags.append("--exclude-atypical")

    if exclude_multi_isolate:
        flags.append("--exclude-multi-isolate")

    if mag != None:
        if mag:
            flags += ["--mag", "only"]
        elif not mag:
            flags += ["--mag", "exclude"]

    return flags


###FUNCIÓN PARA DEFINIR EL ORIGEN DE LOS DATOS
def origen_datos(origen, # taxon | accession
                 nombres #Lista de los valores
                 ):
    if origen not in ["taxon", "accession"]:
        raise ValueError(f"No se ha introducido una de las posibilidades: taxon | accession ({origen})")
    
    ## * permite desempaquetar, para no generar una lista anidada
    return [origen, *nombres]

def descargar_genomas(origen, nombres, filename="ncbi_dataset.zip", flags=None):
    if flags is None:
        flags = []
    
    comando = ["datasets", "download", "genome", *origen_datos(origen,nombres), 
               "--dehydrated", "--filename", filename, *flags]
    
    print(f"Ejecutando: \n {' '.join(comando)}")
    subprocess.run(cmd, check=True)
    return


flags = definir_flags(assembly_source="RefSeq", assembly_level=["complete", "chromosome"],
                       exclude_atypical=True, mag=False)

descargar_genomas("taxon", ["Klebsiella pneumoniae"], filename="kpn.zip", flags=flags)

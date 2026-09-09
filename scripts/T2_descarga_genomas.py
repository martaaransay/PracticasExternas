# Descarga automatizada de genomas con criterios modificables

import subprocess # subprocess permite ejecutar comandos en la terminal (ej. datasets)
import json
import os
#------------------------ Función para definir las flags para la descarga de genomas
def definir_flags(
          assembly_source = None, # RefSeq | GenBank | None =todas
          assembly_level = None, # comma-separated: chromosome | complete | contig | scaffold | None =todas
          annotated = None, # Si no es None, excluye genomas anotados
          exclude_atypical = None, # Si no es None, excluye assemblies atípicos
          exclude_multi_isolate = None, # Si no es None, excluye assemblies de proyectos multi-isolate
          mag = None, # Si es True o 'only', solo usa mags; si es false o 'exclude', los excluye
          ):
    
    # Lista de las flags a rellenar
    flags = []

    if assembly_source: # Si se ha introducido un argumento válido, se añade a la lista de flags
        assembly_source = assembly_source.lower() # Control de mayúsculas/minúsculas
        if assembly_source == "refseq":
            flags += ["--assembly-source", "RefSeq"]
        elif assembly_source == "genbank":
            flags += ["--assembly-source", "GenBank"]
        else:
            raise ValueError(f"No se ha introducido una de las posibilidades: RefSeq | GenBank ({assembly_source})")
        
    if assembly_level: ##FIXME
        ass_level = []
        if isinstance(assembly_level, str):
            assembly_level = [assembly_level]
        for a in assembly_level:
            ass_level.append(a.lower())
        for a2 in ass_level:
            if a2 not in ["chromosome","complete","contig","scaffold"]:
                raise ValueError(f"No se ha introducido una de las posibilidades: chromosome | complete \
                                  | contig | scaffold ({assembly_level})")
        else:
            # Comma-separated:
            flags += ["--assembly-level", ",".join(ass_level)]

    if annotated: # Si se ha introducido el argumento, se añade a la lista de flags
        flags.append("--annotated")

    if exclude_atypical:
        flags.append("--exclude-atypical")

    if exclude_multi_isolate:
        flags.append("--exclude-multi-isolate")

    if mag is not None: # Si se ha introducido un argumento válido, se añade a la lista de flags
        if not isinstance(mag, bool):
            mag = mag.lower()
        if mag is True or mag == "only":
            flags += ["--mag", "only"]
        elif mag is False or mag == "exclude":
            flags += ["--mag", "exclude"]
        elif mag == "all":
            pass
        else:
            raise ValueError(f"No se ha introducido una de las posibilidades: only | exclude | all ({mag})")
    return flags

#------------------------ Función para definir taxon / accession
def origen_datos(origen, # taxon | accession | inputfile
                 nombre # Lista de los valores  o nombre del archivo
                 ):
    origen = origen.lower()
    if origen not in ["taxon", "accession", "inputfile"]:
        raise ValueError(f"No se ha introducido una de las posibilidades: taxon | accession | inputfile ({origen})")
    
    if origen == "inputfile":
        lista = ["accession", "--inputfile", nombre]
    else:
        if isinstance(nombre, str):
            nombre = [nombre] # Para que no separe por letras al desempaquetar
        lista = [origen, *nombre] ## * permite desempaquetar, para no generar una lista anidada
    return lista 


def filtrado_genomas(origen, nombre,
####los umbrales son random ns 
                    umbral_n50=5000,
                    max_contaminacion=2,
                    umbral_continuidad=97):
    # Filtrado de genomas (N50, contaminación...)

    # JSON con los datos:
    comando_summary = ["datasets", "summary", "genome", 
                       *origen_datos(origen, nombre), "--as-json-lines"]

    try:
        summary_json = subprocess.run(comando_summary, 
                                      capture_output = True, # Captura stdout y stderr
                                      text = True, # Los ficheros de stdout y stderr se abren en modo texto
                                      check = True) ### Si check es verdadero y el proceso retorna un resultado distinto de cero, 
                                                    ### se lanzará una excepción CalledProcessError. 
                                                    ### Los atributos de dicha excepción contendrán los argumentos, el código retornado 
                                                    ### y tanto stdout como stderr si se capturaron
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar datasets summary:\n{e.stderr}")
        return

    acc_validos = []
    summary_text = summary_json.stdout.splitlines()
    # print(summary_text)
    for line in summary_text:
        if not line.strip():
            continue
        datos = json.loads(line)
        acc = datos.get("accession")
        # Métricas de calidad:
        #-------------------N50
        assembly_stats = datos.get("assembly_stats", {})
        n50 = assembly_stats.get("scaffold_n50",0) ##FIXME HE PUESTO SCAFFOLD en lugar de contig_n50 pero no se cual seria mejor
        
        #-------------------CONTAMINACIÓN Y CONTINUIDAD
        checkm_info = datos.get("checkm_info", {})
        # Contaminación
        contaminacion = checkm_info.get("contamination", 0.0) # Si no hay datos, asumimos 0 o lo que decidas
        # Continuidad
        continuidad = checkm_info.get("completeness", 100.0)

        
        # Aplicar los filtros
        if n50 >= umbral_n50 and contaminacion <= max_contaminacion and continuidad >= umbral_continuidad:
            print(acc)
            acc_validos.append(acc)

    if not acc_validos:
        print("Ningún genoma cumple con los criterios de calidad.")
        return
        # Guardar los accessions en un archivo de texto temporal
    
    ## FIXME -> UNA OPCIÓN VIABLE SERÍA CAMBIAR EL ARCHIVO POR UNA LISTA TAL CUAL peero HAY UN LÍMITE MÁXIMO DE CARACTERES PARA UN COMANDO EN LA TERMINAL (ARG_MAX)
    archivo_lista = f"lista_accessions_{nombre[0].replace(' ','_')}.txt" #por si acaso hay espacios
    with open(archivo_lista, "w") as f:
        for acc in acc_validos:
            f.write(f"{acc}\n")
    return "inputfile", archivo_lista

##------------------------ Función para descargar los genomas
def descargar_genomas(origen, nombre, 
                      filtrado_calidad = None, # Flag para saber si filtrar por calidad
                      filename="ncbi_dataset.zip", #puede ser un path para especificar el directorio destino
                      flags=None):
    if filtrado_calidad:
        resultado_filtrado = filtrado_genomas(origen, nombre) ##FIXME falta añadir los nuevos argumentos
        if not resultado_filtrado:
            return
        origen, nombre = resultado_filtrado
    # Si no se han añadido flags, lista vacía
    if flags is None:
        flags = []
    ########### -------------------- download
    # Comando completo
    comando_download = ["datasets", "download", "genome", *origen_datos(origen, nombre), 
                        "--dehydrated", "--filename", filename, *flags]
    
    print(f"Ejecutando: \n {' '.join(comando_download)}")
    try:
        subprocess.run(comando_download, 
                        check=True, capture_output=True, text=True) 
    # Capturar la excepción
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar datasets download:\n{e.stderr}")


def descomprimir_rehidratar(filename="ncbi_dataset.zip"): #puede ser un path
    
    if not os.path.isfile(filename):
        print(f"No se ha encontrado el archivo {filename}")
        return
    ############# ------------------- descomprimir

    #se pueden añadir -o (sobreescribe archivos sin preguntar) -q(modo silencioso) / -d (para especificar el directorio de destino)
    # Si filename es solo "archivo.zip", devuelve "" y usamos "." (directorio actual)
    dir_out = os.path.dirname(filename) or "."
    comando_unzip = ["unzip", filename, "-q", "-o", "-d", dir_out]
    print(f"Ejecutando: \n {' '.join(comando_unzip)}") 
    try:
        subprocess.run(comando_unzip, 
                       check=True, capture_output=True, text=True)
    # Capturar la excepción
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar unzip:\n{e.stderr}")
        return
    ############ ------------------- rehydrate

    comando_rehidratar = ["datasets", "rehydrate", "--directory", dir_out]
    print(f"Ejecutando: \n {' '.join(comando_rehidratar)}") 
    try:
        subprocess.run(comando_rehidratar, 
                        check=True, capture_output=True, text=True) 
    # Capturar la excepción
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar datasets rehydrate:\n{e.stderr}")

    return


# a,b =filtrado_genomas("accession", "GCF_000001405.40")

# ##EJEMPLO uso
# flags = definir_flags(assembly_source="RefSeq", assembly_level=["complete", "chromosome"],
#                        exclude_atypical=True, mag=False)

# descargar_genomas("taxon", ["Klebsiella pneumoniae"], filename="kpn.zip", flags=flags)

# descargar_genomas("accession", "GCF_000001405.40", filtrado_calidad=True)


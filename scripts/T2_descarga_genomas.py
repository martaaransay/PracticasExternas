#### T2: Descarga automatizada de genomas desde NCBI
# Añadir criterios modificables

#### LIBRERÍAS:
import subprocess # para ejecutar comandos en la terminal (ej. datasets / unzip)
import shlex # para comprobar qué se envía exactamente a la consola
import json # para poder utilizar correctamente .json
import os # comprobación de existencia de ciertos archivos

#### ------- FUNCIONES:

#------------------------ Función para definir las flags para la descarga de genomas (datasets download)
def definir_flags(
          assembly_source = None, # RefSeq | GenBank | None =todas
          assembly_level = None, # comma-separated (string): chromosome | complete | contig | scaffold | None =todas
          annotated = None, # Para excluir genomas anotados (≠ None)
          exclude_atypical = None, # Excluye assemblies atípicos (≠ None) #FIXME -> comprobar qué considera como atípicos
          exclude_multi_isolate = None, # Excluye assemblies de proyectos multi-isolate (≠ None) #FIXME
          mag = None, # Si es True o 'only', solo usa MAGS; si es false o 'exclude', los excluye
          ):
    """
    Construye y valida una lista de flags para el comando de descarga de NCBI datasets.

    Args
    --------------------
    assembly_source (str): Origen del ensamblaje. Valores permitidos: 'RefSeq' o 'GenBank'.
    assembly_level (str | list): Nivel de ensamblaje. Valores permitidos: 'chromosome', 'complete', 'contig', 'scaffold'.
    annotated (bool): Si evalúa a ≠ False/None, incluye genomas anotados ('--annotated').
    exclude_atypical (bool): Si evalúa a ≠ False/None, excluye ensamblajes atípicos ('--exclude-atypical').
    exclude_multi_isolate (bool): Si evalúa a ≠ False/None, excluye proyectos multi-aislado ('--exclude-multi-isolate').
    mag (bool | str): Filtro para genomas metagenómicos. Valores: True/'only', False/'exclude' o 'all'.

    Returns
    ----------------------
    list: Lista de cadenas de texto con los argumentos listos para ser pasados a subprocess.
"""
    
    flags = [] # Lista de las flags a rellenar

    if assembly_source: # Si se ha introducido un argumento válido, se añade a la lista de flags
        assembly_source = assembly_source.lower() # Control de mayúsculas/minúsculas
        if assembly_source == "refseq":
            flags += ["--assembly-source", "RefSeq"] # Se añade a la lista el argumento 
        elif assembly_source == "genbank":
            flags += ["--assembly-source", "GenBank"]
        else:
            raise ValueError(f"No se ha introducido una de las posibilidades: RefSeq | GenBank \
                            \n (introducido: {assembly_source})")
        
    if assembly_level: 
        ass_level = []
        if isinstance(assembly_level, str): # Por si se ha introducido un string:
            for item_str in assembly_level.split(","): #comma-separated
                ass_level_indiv = item_str.strip().lower() # Control de blancos y mayúsculas
                ass_level.append(ass_level_indiv) 
        elif isinstance(assembly_level, list): # Por si se ha introducido una lista:
            for item_list in assembly_level:
                ass_level_indiv = item_list.strip().lower()
                ass_level.append(ass_level_indiv)
        else: # Formato incorrecto
            raise ValueError(f"No se ha introducido el formato correcto (lista o string comma-separated) \
                                    \n (introducido: {assembly_level})")
        
        for item in ass_level: # Control de argumentos válidos
            if item not in ["chromosome","complete","contig","scaffold"]:
                raise ValueError(f"No se ha introducido una de las posibilidades: chromosome | complete | contig | scaffold \
                                \n (introducido: {assembly_level})")
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
        if not isinstance(mag, bool): # Control de mayúsculas para no-booleanos
            mag = mag.lower()
        if mag is True or mag == "only":
            flags += ["--mag", "only"]
        elif mag is False or mag == "exclude":
            flags += ["--mag", "exclude"]
        elif mag == "all":
            pass
        else:
            raise ValueError(f"No se ha introducido una de las posibilidades: only | exclude | all \
            \n (introducido: {mag})")
    return flags


#------------------------ Función para definir taxon / accession
def origen_datos(origen, # Args válidos: taxon | accession | inputfile
                 nombre # Lista de los valores (origen="taxon" | origen="accession")\
                        # o nombre del archivo (origen="inputfile")
                 ):
    """
    Construye la lista de argumentos respecto al los genomas a descargar para el comando de descarga de NCBI datasets.

    Args
    -----------------
    origen (str): Tipo de entrada para la búsqueda. Valores permitidos: 'taxon', 'accession' o 'inputfile'.
    nombre (str | list): Valores a buscar (string separado por comas o lista de strings) o ruta del archivo si el origen es 'inputfile'.

    Returns
    -----------------
    list: Lista de cadenas de texto con los argumentos del origen desempaquetados, listos para el comando de subprocess.
    """
    origen = origen.lower() # Control de mayúsculas
    if origen not in ["taxon", "accession", "inputfile"]: # Control de argumentos válidos
        raise ValueError(f"No se ha introducido una de las posibilidades: taxon | accession | inputfile \
                        \n (introducido: {origen})")
    
    if origen == "inputfile":# File
        if not os.path.isfile(nombre): # Control por si no existe el archivo
            print(f"No se ha encontrado el archivo {nombre}")
            return
        lista = ["accession", "--inputfile", nombre] # Nombre del archivo directamente
    
    else:
        nombre_nofile=[]
        if isinstance(nombre, str): # Elemento único o varios comma-separated
            for item_nombre_str in nombre.split(","): # Si es elemento único (no ,) no lanza error
                # El control de mayúsculas no es necesario (probado en datasets)
                nombre_nofile.append(item_nombre_str)
        elif isinstance(nombre, list): # Lista
            for item_nombre_lista in nombre:
                nombre_nofile.append(item_nombre_lista)
        else:# Formato incorrecto
            raise ValueError(f"No se ha introducido el formato correcto (lista o string comma-separated) \
                                    \n (introducido: {nombre})")

        lista = [origen, *nombre_nofile] ## * permite desempaquetar, para no generar una lista anidada
    return lista 


def filtrado_genomas(origen, nombre, # Argumentos para la función origen_datos(origen, nombre)
                    filename, # Filename para mantener el mismo directorio
                    # FIXME -> Seleccionar umbrales válidos 
                    umbral_n50=5000, # Umbral para el N50
                    max_contaminacion=2, # Umbral de contaminación
                    umbral_continuidad=97): # Umbral de continuidad 
    """
    Consulta los metadatos de genomas en NCBI datasets, filtra por métricas de calidad y guarda los accessions válidos en un fichero de texto.

    Args
    ------------
    origen (str): Tipo de origen para la consulta ('taxon', 'accession' o 'inputfile').
    nombre (str | list): Valor o identificador asociado al origen (taxón, lista de accessions o ruta de fichero).
    filename (str): Nombre del archivo .zip resultante. 
    umbral_n50 (int, opcional): Valor mínimo aceptable para el scaffold N50. Por defecto es 5000.
    max_contaminacion (float | int, opcional): Porcentaje máximo tolerado de contaminación según CheckM. Por defecto es 2.
    umbral_continuidad (float | int, opcional): Porcentaje mínimo requerido de completitud (completeness) según CheckM. Por defecto es 97.

    Return
    --------------
    tuple[str, str] | None: Tupla ('inputfile', ruta_del_archivo) con el nombre del fichero generado si hay genomas que superan los filtros
                            None si no hay resultados válidos o si la ejecución falla.
    """
    # Para obtener el JSON con los datos sobre el genoma 
    comando_summary = ["datasets", "summary", "genome", 
                       *origen_datos(origen, nombre), "--as-json-lines"]
    print(f"Ejecutando: \n {shlex.join(comando_summary)}")
    try:
        summary_json = subprocess.run(comando_summary, 
                                      capture_output = True, # Captura stdout y stderr
                                      text = True, # ficheros de stdout y stderr en modo texto
                                      check = True) # Para capturar si ocurre un error (excepción CalledProcessError)
    except subprocess.CalledProcessError as e: # Captura la excepción
        print(f"Error al ejecutar datasets summary:\n{e.stderr}")
        return

    acc_validos = [] # Lista para almacenar los accessions válidos 
    summary_text = summary_json.stdout.splitlines() # Json obtenido (separado por líneas)
    for line in summary_text:
        if not line.strip(): continue # Salta líneas vacías
        datos = json.loads(line) # Deserialización
        acc = datos.get("accession") # Obtener el accession
        
        # Métricas de calidad:
        #-------------------N50
        assembly_stats = datos.get("assembly_stats", {}) # Obtener stats de forma segura (si no están, dic vacío)
        n50 = assembly_stats.get("scaffold_n50",0) # FIXME (scaffold_n50 o contig_n50)
        
        #-------------------CONTAMINACIÓN Y CONTINUIDAD
        checkm_info = datos.get("checkm_info", {})
        # Contaminación
        contaminacion = checkm_info.get("contamination", 0.0) # FIXME -> Si no hay información, se asume que no hay contaminación
        # Continuidad
        continuidad = checkm_info.get("completeness", 100.0) # FIXME -> Si no hay información, se asume continuidad 100

        # Comprobar que cumpla los filtros seleccionados 
        if n50 >= umbral_n50 and contaminacion <= max_contaminacion and continuidad >= umbral_continuidad:
            acc_validos.append(acc)

    if not acc_validos:
        print("Ningún genoma cumple con los criterios de calidad.")
        return
    
    ## FIXME -> UNA OPCIÓN VIABLE para no crear un archivo:
        # ES DECIR: CAMBIAR EL ARCHIVO POR UNA LISTA TAL CUAL peero HAY UN LÍMITE MÁXIMO DE CARACTERES PARA UN COMANDO EN LA TERMINAL (ARG_MAX)
    archivo_lista = f"lista_accessions_{os.path.basename(filename).split('.')[0]}.txt"
    ruta_archivo = os.path.dirname(filename) + "/" + archivo_lista
    with open(ruta_archivo, "w") as f: # Creación de un archivo con los acc válidos (uno por línea)
        for acc in acc_validos:
            f.write(f"{acc}\n")
        print(f"Se han guardado en el archivo '{ruta_archivo}' los accessions que se van a analizar.")
    return "inputfile", ruta_archivo


##------------------------ Función para descargar los genomas
def descargar_genomas(origen, nombre, # Argumentos para las funciones: origen_datos() | filtrado_genomas()
                      flag_filtrado_calidad = None, # Flag para saber si filtrar por métricas de calidad
                      filename="ncbi_dataset.zip", #Nombre del archivo que va a tener los genomas 
                                # Puede tener un path para especificar un directorio destino
                      flags = None): #flags para el comando datasets
    """
    Ejecuta el comando de NCBI datasets para descargar genomas (--dehydrated), permitiendo un filtrado previo de calidad.

    Args
    --------
    origen (str): Tipo de entrada original para la descarga ('taxon', 'accession' o 'inputfile').
    nombre (str | list): Identificadores de los genomas, nombre del taxón o ruta del archivo base.
    flag_filtrado_calidad (bool, opcional): Si evalúa a True, ejecuta un filtrado previo de calidad; descartando los que no lo superen.
    filename (str, opcional): Nombre y ruta del archivo .zip resultante. 
                              Por defecto es "ncbi_dataset.zip".
    flags (list, opcional): Lista de flags adicionales para añadir al comando 'datasets download'. 
                            Por defecto es None.
    """
    # Control de que los directorios existan para evitar errores
    dir_out = os.path.dirname(filename)
    if dir_out: # Evita errores si la ruta no tiene carpetas 
        os.makedirs(dir_out, exist_ok=True) 

    if flag_filtrado_calidad:
        resultado_filtrado = filtrado_genomas(origen, nombre, filename) # Ejecuta la función para filtrar
        if not resultado_filtrado: # Si no hay accession, no se puede utilizar el comando datasets
            # Se eliminan los directorios que se han creado:
            if dir_out and os.path.exists(dir_out):
                if not os.listdir(dir_out): 
                    os.removedirs(dir_out) # Elimina directorios vacíos (parentales también)
            return
        origen, nombre = resultado_filtrado # Desempaquetar tupla
    # Si no se han añadido flags, lista vacía (control)
    if flags is None:
        flags = []

    ### Comandos:
    ########### -------------------- datasets download
    comando_download = ["datasets", "download", "genome", *origen_datos(origen, nombre), 
                        "--dehydrated", "--filename", filename, *flags]
    
    print(f"Ejecutando: \n {shlex.join(comando_download)}")
    try:
        subprocess.run(comando_download, 
                       check=True, capture_output=True, text=True) 
    # Capturar la excepción
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar datasets download:\n{e.stderr}")
        # Se eliminan los directorios que se han creado (siempre que estén completamente vacíos)
        if dir_out and os.path.exists(dir_out):
            if not os.listdir(dir_out): 
                os.removedirs(dir_out) # Elimina directorios vacíos (parentales también)
        return None
    return True

##------------------------ Función para descomprimir y rehydrate los genomas
def unzip_rehydrate(filename = "ncbi_dataset.zip"): 
    
    if not os.path.isfile(filename): # Control de que el archivo exista
        print(f"No se ha encontrado el archivo {filename}")
        return

    ############# ------------------- unzip
    dir_out = os.path.dirname(filename) or "." # Si no hay path en el filename, devuelve "."
    nueva_carpeta = os.path.basename(filename)
    dir_out = dir_out + "/" + nueva_carpeta.split(".")[0]
    comando_unzip = ["unzip",
                     "-q",  # Quiet
                     "-o",  # Sobreescribe archivos para evitar errores
                     filename,
                     "-d", dir_out] # Directorio donde se descomprime
    print(f"Ejecutando: \n {shlex.join(comando_unzip)}")
    try:
        subprocess.run(comando_unzip, 
                       check=True, capture_output=True, text=True)
    # Capturar la excepción
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar unzip:\n{e.stderr}")
        return

    ############ ------------------- datasets rehydrate
    comando_rehydrate= ["datasets", "rehydrate", "--directory", dir_out]
    print(f"Ejecutando: \n {shlex.join(comando_rehydrate)}")
    try:
        subprocess.run(comando_rehydrate, 
                        check=True, capture_output=True, text=True) 
    # Capturar la excepción
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar datasets rehydrate:\n{e.stderr}")
    return

####### EJEMPLOS FINALES PROBADOS:

flags = definir_flags(assembly_source="RefSeq", 
                      assembly_level=["complete", "chromosome"],
                      exclude_atypical=True, 
                      mag=False)

file = "results/pruebas_T2/kpn.zip"
check = descargar_genomas(origen="taxon", 
                          nombre=["Klebsiella pneumoniae"], 
                          filename=file, 
                          flags=flags)

if check:
    unzip_rehydrate(filename=file)
# FIXME -> se puede cambiar facilmente para en vez de ejecutar decargar y unzip+ rehydrate hacerlo junto, pero bueno

# FIXME -> no salen las barras de progreso al ejecutar este script, pero si quitamos el argumento (capture_output=True)
#          en subprocess.run, sale -> solo que ya no se captura el error así (habría que eliminar esa parte de:
#          (except subprocess.CalledProcessError as e:)
#          PERO OJO -> SOLO EN LAS QUE HAGA FALTA (EJ. en unzip no!)
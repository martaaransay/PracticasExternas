import os
import subprocess 
import shlex 

def adaptar_dir_IntegronFinder(acc_genoma, directorio):
    """
    Busca y extrae las rutas de los archivos gbk generados por IntegronFinder2.

    Args
    --------
    acc_genoma (str): Nombre de accesión del genoma (para identificar la carpeta). 
    directorio (str): Ruta al directorio donde se encuentran los resultados.

    Returns
    --------
    archivos_gbk (list): Lista que contiene las rutas completas a los archivos .gbk
    """
    archivos_gbk = [] # Lista para almacenar las rutas a los archivos .gbk
    # Identifica el directorio concreto donde se encuentran los resultados para ese genoma
    path_resultados = os.path.join(directorio, f"Results_Integron_Finder_{acc_genoma}")
    resultados_integronfinder = os.listdir(path_resultados) # Lista los archivos en la carpeta
    for archivo_IntegronFinder in resultados_integronfinder:
        if archivo_IntegronFinder.endswith(".gbk"): # Si es un archivo gbk, se almacena su ruta en la lista
            gbk = os.path.join(path_resultados, archivo_IntegronFinder)
            archivos_gbk.append(gbk)
    return archivos_gbk
            
def Integron_Finder(ruta_genoma, # Ruta al genoma
                    acc_genoma, # Nombre del archivo del genoma
                    outdir,  # Directorio de salida de los resultados
                    # Argumentos booleanos:
                    localmax = False, # Argumento de búsqueda más sensible (--localmax)
                    palindromes = False, # Mantener las versiones palindrómicas (--keep-palindromes)
                    eagle_eyes = False, # Argumento de búsqueda más sensible (--eagle_eyes)
                    # Argumentos con valor path:
                    attc_model = None, # Otro modelo de attC a usar (--attc-model ATTC_MODEL)
                    topology_file = None, # Ruta al archivo con las topologías de los replicones [--topology-file TOPOLOGY_FILE]
                    # Argumentos con valor numérico:
                    cpu = min(os.cpu_count()-1, 4), # Número de CPUs a utilizar (máx 4 por recomendación de integronfinder2) (--cpu)
                    distance_threshold = 4000, # Threshold de distancia, default en IntegronFinder en 4000 (--dt DISTANCE_THRESHOLD)
                    evalue_attc = 1, # E-value usado para filtrar attC (--evalue-attc EVALUE_ATTC)
                    max_attc_size = 200, # Máximo tamaño del sitio attC (--max-attc-size MAX_ATTC_SIZE)
                    min_attc_size = 40, # Mínimo tamaño del sitio attC (--min-attc-size MIN_ATTC_SIZE)
                    ): 
    """
    Ejecuta IntegronFinder2 sobre un genoma de entrada con parámetros personalizables.

    Args
    -------
    ruta_genoma (str): Ruta completa al archivo del genoma de entrada.
    acc_genoma (str): Nombre de accesión del genoma 
    outdir (str): Ruta al directorio donde se guardarán los resultados.
    localmax (bool, opcional): Activa una búsqueda más sensible de sitios attC (--local-max). 
                               Por defecto False.
    palindromes (bool, opcional): Mantiene las versiones palindrómicas de los sitios attC (--keep-palindromes). 
                                  Por defecto False.
    eagle_eyes (bool, opcional): Reduce el umbral de ruido para buscar elementos difíciles de detectar (--eagle-eyes). 
                                 Por defecto False.
    attc_model (str, opcional): Ruta a un modelo de covarianza alternativo para los sitios attC (--attc-model).
                                Por defecto None.
    topology_file (str, opcional): Ruta al archivo que define las topologías (circular/lineal) de los replicones (--topology-file). 
                                   Por defecto None.
    cpu (int, opcional): Número de núcleos de CPU a utilizar. 
                         Por defecto se usan 4 núcleos como máximo.
    distance_threshold (int, opcional): Distancia máxima en pares de bases entre dos elementos de un integrón (-dt). 
                                        Por defecto 4000.
    evalue_attc (float o int, opcional): Umbral de E-value para filtrar los sitios attC encontrados (--evalue-attc). 
                                         Por defecto 1.
    max_attc_size (int, opcional): Tamaño máximo permitido para un sitio attC (--max-attc-size). 
                                   Por defecto 200.
    min_attc_size (int, opcional): Tamaño mínimo permitido para un sitio attC (--min-attc-size). 
                                   Por defecto 40.

    Return
    -------
    archivos_gbk (list): Lista de los archivos .gbk resultantes tras ejecutar IntegronFinder2.
                         Devuelve None si ocurre un error en la ejecución.
    """

    ### Comando general IntegronFinder2
    comando_integronfinder = ["integron_finder", ruta_genoma,
                              "--promoter-attI", # --promoter-attI para buscar promotores
                              "--outdir", outdir, "--gbk", # --gbk para el output
                              "--cpu", str(cpu), "-dt", str(distance_threshold),
                              "--evalue-attc", str(evalue_attc), 
                              "--max-attc-size", str(max_attc_size),
                              "--min-attc-size", str(min_attc_size)] 
    
    # Adición de argumentos booleanos:
    if localmax:
        comando_integronfinder.append("--local-max")
    if palindromes:
        comando_integronfinder.append("--keep-palindromes")
    if eagle_eyes:
        comando_integronfinder.append("--eagle-eyes")

    # Adición de argumentos con valor path:
    if attc_model is not None:
        if not os.path.isfile(attc_model): # Control de path correcto
            print(f"attc_model no encontrado: {attc_model}")
        else:
            comando_integronfinder += ["--attc-model", attc_model]
        
    if topology_file is not None:
        if not os.path.isfile(topology_file): # Control de path correcto
            print(f"topology_file no encontrado: {topology_file}")
        else:
            comando_integronfinder += ["--topology-file", topology_file]

    print(f"Ejecutando: \n {shlex.join(comando_integronfinder)}")
    try:
        subprocess.run(comando_integronfinder, 
                       check=True, capture_output=True, text=True)
    # Capturar la excepción
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar IntegronFinder2 :\n{e.stderr}")
        return None
        
    archivos_gbk = adaptar_dir_IntegronFinder(acc_genoma, outdir)
    return archivos_gbk


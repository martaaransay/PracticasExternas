import os
from T3_0_integronfiltering import integron_filtering
from T3_1_integronfinder import Integron_Finder
from T3_2_Pc import extraer_Pc_gbk

def ajustar_directorio_T2(dir_T2):
    """
    Adapta el directorio indicado en el anterior script al necesario para continuar el pipeline.
    
    Args
    -------
    dir_T2 (str): Directorio dado en el script anterior

    Return
    ------
    dir_adaptado (str): Nuevo directorio adaptado
    """
    archivo = os.path.basename(dir_T2) # Extrae el nombre y extensión del archivo
    nombre_archivo = os.path.splitext(archivo)[0] # Extrae únicamente el nombre del archivo
    # Añade al nuevo directorio el nombre del archivo del genoma y las carpetas generadas al utilizar datasets
    dir_adaptado = os.path.join(os.path.dirname(dir_T2), 
                               nombre_archivo, 
                               "ncbi_dataset/data") # Estructura generada por datasets
    return dir_adaptado

def acceder_archivos_fasta(directorio):
    """
    Obtiene las rutas de todos los archivos FASTA (.fna) almacenados dentro de las subcarpetas de un directorio.

    Args
    -------
    directorio (str): Ruta del directorio principal que contiene las subcarpetas de los genomas.

    Returns
    --------
    rutas_genomas (list): Una lista con las rutas completas hacia cada archivo fasta (.fna) encontrado.
    """
    rutas_genomas = [] # Lista para almacenar la ruta de cada genoma
    carpetas = sorted(os.listdir(directorio)) # lista de carpetas dentro del directorio adaptado
    for carpeta in carpetas: 
        ruta_individual = os.path.join(directorio, carpeta) # Añade la carpeta dentro de la ruta
        # Si se ha descargado correctamente, se genera una nueva carpeta que contiene dentro el genoma
        if os.path.isdir(ruta_individual): 
            archivos_por_carpeta = os.listdir(ruta_individual) # Lista de los archivos de la nueva carpeta
            for archivo in archivos_por_carpeta:
                if os.path.splitext(archivo)[-1] == ".fna": 
                    # Si el archivo concreto se corresponde con un .fasta, se añade a la ruta
                    ruta_archivo = ruta_individual + "/" + archivo
                    rutas_genomas.append(ruta_archivo) # Guarda la ruta individual de cada genoma
    return rutas_genomas

def recorrer_genomas(rutas_genomas,
                     directorio,
                     flag_integronfiltering = True): # Por defecto se filtra previamente con integron_filtering
    """
    FIXME -> Falta docstring cuando se acabe la función
    """
    dir_parental = directorio.replace("/ncbi_dataset/data", "") # Utiliza como referencia el directorio parental 
    
    for ruta_genoma_individual in rutas_genomas: # Recorre cada archivo con los genomas
        archivo = os.path.basename(ruta_genoma_individual) # Extrae el nombre y extensión del archivo
        acc_genoma = os.path.splitext(os.path.basename(archivo))[0] # Extrae únicamente el nombre del archivo
        # Filtrado vía integron-filtering
        if flag_integronfiltering:
            # Crea una carpeta para los resultados del filtrado, y dentro de esta un subdirectorio con cada genoma         
            directorio_integronfiltering = os.path.join(dir_parental, "results_IntegronFiltering", acc_genoma)
            # Crea el directorio (no error si ya existe)
            os.makedirs(directorio_integronfiltering, exist_ok=True) 
            # Filtrado vía integron-filtering y actualización de la ruta con el genoma filtrado
            ruta_genoma_individual = integron_filtering(ruta_genoma_individual, directorio_integronfiltering)
            # Actualización del accession filtrado
            archivo = os.path.basename(ruta_genoma_individual)
            acc_genoma = os.path.splitext(os.path.basename(archivo))[0]
            
        # Crea una carpeta para los resultados de la búsqueda      
        directorio_integronfinder = os.path.join(dir_parental, "results_IntegronFinder2")
        # Crea el directorio (no error si ya existe)
        os.makedirs(directorio_integronfinder, exist_ok=True) 
        # Ejecución de IntegronFinder2
        archivos_gbk_integronfinder = Integron_Finder(ruta_genoma_individual, 
                                                     acc_genoma, 
                                                     directorio_integronfinder)
        if not archivos_gbk_integronfinder:
            print(f"IntegronFinder2 no ha detectado integrones en el genoma {ruta_genoma_individual}.")
            # FIXME: Habría que borrar aquí algún archivo / directorio para que no de problemas al extraer PCs y clasificar??
            continue
            
        ## EXTRAER INFORMACIÓN SOBRE PCs
        lista_pc = extraer_Pc_gbk(acc_genoma, archivos_gbk_integronfinder)
        for element in lista_pc:
            print(element["secuencia"])
        
        ### falta generar el módulo de clasificación de Pcs!


### FIXME -> Falta arreglar la parte de que no salga Promotor en los .gbk por ejemplo!!!!
### FIXME -> Muchas veces aparece además el mensaje de: 
# "No se ha podido procesar el genoma results/pruebas_T2/kpn/results_IntegronFiltering/GCF_001663295.1_ASM166329v1_genomic/GCF_001663295.1_ASM166329v1_genomic.attC_filtered.fasta por IntegronFinder2"


# EJEMPLOS USADOS:


path=ajustar_directorio_T2("results/pruebas_T2/kpn.zip")

genom = acceder_archivos_fasta(path)


recorrer_genomas(genom, path)

# hacer pruebas con flag de filtrado False
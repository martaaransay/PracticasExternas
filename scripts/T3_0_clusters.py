import subprocess
import shlex
import os
import glob # para interpretar correctamente los *

### LUEGO BORRAR ESTA FUNCION
file = "results/T2_check/kpn_filtered.zip"
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

def cluster_genomes(genomes_dir,
                    out_dir,
                    threshold = 0.99):
    """docstring"""

    # Check if the directory exists

    if not os.path.isdir(genomes_dir):
        raise FileNotFoundError(f"The directory {genomes_dir} was not found.")
    
    genomes_path = os.path.join(genomes_dir, "*", f"*.fna")
    genomes_path_glob = glob.glob(genomes_path)

    dRep_cmd = ["dRep", "dereplicate", out_dir,
                "-g", *genomes_path_glob, "-sa", str(threshold), 
                #FIXME: He tenido que añadir estos argumentos porque sino petaba, pporbar en el otro ordenador
                "--multiround_primary_clustering", "--primary_chunksize", "500", 
                "--low_ram_primary_clustering", "-p", "4",
                "--ignoreGenomeQuality"] 
    # print(f"Executing: \n{shlex.join(dRep_cmd)}")
    try:
        a = subprocess.run(dRep_cmd, 
                                      capture_output = True, # Captures stdout and stderr
                                      text = True, # stdout and stderr files in text mode
                                      check = True) # To capture if an error occurs (CalledProcessError exception)
    except subprocess.CalledProcessError as e: # Captures the exception
        print(f"Error executing datasets summary:\n{e.stderr}")
        return

file_general = "results/T2_check/Deinococcus.zip"
file_adaptado = ajustar_directorio_T2(file_general)
cluster_genomes(file_adaptado, "results/T2_check/clustering99")

## FIXME: cambiar lo de -g a: Can also input a text file with paths to genomes, which results in fewer OS issues than wildcard expansion (default: None)
# Import libraries 
import subprocess
import shlex
import os
import glob # para interpretar correctamente los *
import csv

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

# --------Definition of functions
def run_checkm2(genomes_path_glob, out_dir, 
                lowmem = None, 
                threads = 4):
    """
    docstring
    """
    #FIXME: poner en el readme o así del script la creación del entorno conda:
    # subprocess.run(["conda", "create", "-n", "checkm2", "-c", "conda-forge", "checkm2", "-y"], check=True)
    ## FIXME: Y TMB LA DATABASE DE DIAMOND
    # subprocess.run(["conda", "run", "-n", "checkm2", "checkm2", "database", "--download"], check=True)
    
    checkm2_cmd = ["checkm2", "predict", "--threads", str(threads),
                   "--input", *genomes_path_glob, "-x", "fna",
                   "--output-directory", out_dir,
                   "--remove_intermediates"]

    if lowmem: # If an argument was provided, adds it to the flags list
        checkm2_cmd.append("--lowmem")

    
    comando_final = ["conda", "run", "-n", "checkm2", *checkm2_cmd]
 
    try:
        subprocess.run(comando_final, capture_output=True, text=True, check=True)
    
    except subprocess.CalledProcessError as e:
        print(f"Error executing checkm2 predict:\n{e.stderr}")
        return None
    
    return os.path.join(out_dir, "quality_report.tsv")


def checkm2_to_genomeinfo(quality_report_tsv, genome_ext, out_csv):
    """
docstring
    """
    with open(quality_report_tsv) as f_in, open(out_csv, "w", newline="") as f_out:
        reader = csv.DictReader(f_in, delimiter="\t")
        writer = csv.writer(f_out)
        writer.writerow(["genome", "completeness", "contamination"])
        for row in reader:
            genome_name = f"{row['Name']}.{genome_ext}"
            writer.writerow([genome_name, row["Completeness"], row["Contamination"]])
    return out_csv


def cluster_genomes(genomes_dir,
                    out_dir,
                    threshold = 0.99,
                    threads = 4,
                    length = 50000,
                    completeness = 75,
                    contamination = 25,
                    skani_extra = None,
                    cov_thresh = 0.1,
                    chunksize = 500):
    """docstring"""

    if not os.path.isdir(genomes_dir):
        raise FileNotFoundError(f"The directory {genomes_dir} was not found.")

    genomes_path = os.path.join(genomes_dir, "*", "*.fna")
    genomes_path_glob = glob.glob(genomes_path)

    os.makedirs(out_dir, exist_ok=True)
    
    checkm2_out = os.path.join(out_dir, "checkm2")
    quality_report = run_checkm2(genomes_path_glob, checkm2_out, threads=threads)
    print(quality_report)
    if quality_report is None:
        print("CheckM2 falló, abortando clustering.")
        return
    return
    genomeinfo_csv = checkm2_to_genomeinfo(quality_report, "fna",
                                            os.path.join(out_dir, "genomeInfo.csv"))

    # Fichero de texto con rutas en vez de wildcard expandido (recomendado por dRep
    # para evitar problemas de OS con miles de genomas)
    genomes_list_file = os.path.join(out_dir, "genomes_list.txt")
    with open(genomes_list_file, "w") as f:
        f.write("\n".join(genomes_path_glob))
#####FIXME FALTA CAMBIAR TODO ESTO!
    dRep_cmd = ["dRep", "dereplicate", "-p", str(threads),
                "-l", str(length), "-comp", str(completeness), "-con", str(contamination), 
                "--genomeInfo", genomeinfo_csv, "--S_algorithm", "skani", 
                "-sa", str(threshold), "-nc", str(cov_thresh), 
                "--low_ram_primary_clustering", "--multiround_primary_clustering", 
                "-g", genomes_list_file,  "--primary_chunksize", str(chunksize), 
                "--greedy_secondary_clustering", "--run_tertiary_clustering", 
                "--skip_plots",  out_dir
                ]
    
    if skani_extra: # FIXME: No se si está bien, igual con * o con +=
        dRep_cmd.append(skani_extra)

    try:
        subprocess.run(dRep_cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing dRep dereplicate:\n{e.stderr}")
        return


file_general = "results/T2_check/Deinococcus.zip"
file_adaptado = ajustar_directorio_T2(file_general)
cluster_genomes(file_adaptado, "results/T2_check/clustering99")
import subprocess 
import shlex
import os

# FIXME: Cuando acabe, hay que adaptar esto a cómo se encuentre finalmente!
actual_dir = os.path.dirname(__file__) # __file__ ruta del .py actual ()
HattCI_dir = os.path.join(actual_dir, "T3_integronfiltering", "HattCI") # Ruta donde está el ejecutable de HattCI
# os.environ es un dict con las variables de entorno -> os.environ["PATH"] devuelve la variable PATH
# os.pathsep es el caracter separador de PATH en cada sistema operativo
os.environ["PATH"] = f"{HattCI_dir}{os.pathsep}{os.environ['PATH']}" # Sobreescribe la variable PATH añadiendo el ejecutable de HattCI
# Al poner antes el diretorio de HattCI, si había alguno más, se prioriza el indicado antes
#### !!!! el cambio de PATH solo dura mientras se ejecuta el script


# ---------------- FUNCIÓN PARA EJECUTAR INTEGRON FILTERING

def adaptar_dir_integronfiltering(outdir):
    """
    Busca y devuelve la ruta del archivo FASTA filtrado generado por integron_filtering.

    Args
    ------
    outdir (str): Ruta al directorio que contiene los resultados generados por el 
                  script de filtrado de integrones.

    Return:
    ------
    nuevo_archivo (str): La ruta completa al archivo que contiene el genoma con contexto.
    """
    # Archivos generados por el integron_filtering .sh
    resultados_integronfiltering = os.listdir(outdir)
    for resultado_archivo in resultados_integronfiltering:
        # El archivo con el genoma con contexto (NO solo los cortes) se almacena como nueva ruta
        if resultado_archivo.endswith("attC_filtered.fasta"):
            nuevo_archivo = os.path.join(outdir, resultado_archivo)
            return nuevo_archivo


def integron_filtering(ruta_genoma, # Ruta al genoma
                       outdir, # Directorio de salida de los resultados
                       modelo_covarianza = os.path.join(actual_dir, "T3_integronfiltering/attCs.cm"), # default, modelo covarianza definido
                       cutoff_bit = 20, # cut-off del bit score del modelo de covarianza (default en el .sh)
                       min_pb = 500, # mínimo tamaño de la secuencia reconocida (default en el .sh)
                       cpu = os.cpu_count() -1 # número de CPUs a utilizar (cambiado, default en el .sh es 1)
                       ):
    """
    Ejecuta el script de filtrado de integrones (attC-screening).

    Args
    ------
    ruta_genoma (str): Ruta al archivo del genoma que se va a procesar.
    outdir (str): Directorio de salida donde se guardarán los resultados.
    modelo_covarianza (str, opcional): Ruta al archivo del modelo de covarianza. 
                                       Por defecto es 'T3_integronfiltering/attCs.cm'.
    cutoff_bit (int/float, opcional): Valor de corte (cut-off) del bit score para el filtrado. 
                                      Por defecto es 20.
    min_pb (int, opcional): Tamaño mínimo en pares de bases de la secuencia reconocida. 
                            Por defecto es 500.
    cpu (int, opcional): Número de CPUs a utilizar para el proceso. 
                         Por defecto usa el total disponible menos 1.

    Return
    ------
    ruta_genoma_filtrado (str): La ruta al genoma filtrado generado si el proceso es exitoso. 
                                Si ocurre un error durante la ejecución de subprocess, devuelve la ruta original del genoma (ruta_genoma).
    """
    # Comando integron_filtering 
    comando_integronfiltering = ["./scripts/T3_integronfiltering/attC-screening.sh", # ruta al .sh
                                 "-m", modelo_covarianza, 
                                 "-i", ruta_genoma,
                                 "-o", outdir,
                                 # Los números deben pasarse como str
                                 "-b", str(cutoff_bit), 
                                 "-l", str(min_pb),
                                 "-t", str(cpu)]

    print(f"Ejecutando: \n {shlex.join(comando_integronfiltering)}")
    try:
        subprocess.run(comando_integronfiltering, 
                       capture_output = True, text = True, check = True) 
    # Capturar la excepción
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar integron_filtering :\n{e.stderr}")
        return ruta_genoma # Ruta inicial
    # Función para adaptar la ruta al nuevo genoma filtrado
    ruta_genoma_filtrado = adaptar_dir_integronfiltering(outdir)
    # Control por si no se ha generado el archivo correctamente
    if not ruta_genoma_filtrado:
        print(f"Error al ejecutar integron_filtering")
        return ruta_genoma
    
    return ruta_genoma_filtrado

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
#### ojo -> el cambio de PATH solo dura mientras se ejecuta el script


# ---------------- FUNCIÓN PARA EJECUTAR INTEGRON FILTERING

def adaptar_dir_integronfiltering(outdir):
    resultados_integronfiltering = os.listdir(outdir)
    for resultado_archivo in resultados_integronfiltering:
        if resultado_archivo.endswith("attC_filtered.fasta"):
            nuevo_archivo=os.path.join(outdir,resultado_archivo)
            return nuevo_archivo


def integron_filtering(genoma,
                       outdir,
                       modelo_covarianza = "./scripts/T3_integronfiltering/attCs.cm", # default, modelo covarianza definido
                       cutoff_bit = 20, # default en el propio .sh
                       min_pb = 500, # default en el propio .sh
                       cpu=1 # default en el propio .sh
                       ):

    comando_integronfiltering = ["./scripts/T3_integronfiltering/attC-screening.sh", 
                                 "-m", modelo_covarianza,
                                 "-i", genoma,
                                 "-o", outdir,
                                 "-b", str(cutoff_bit), # deben ser str no int!
                                 "-l", str(min_pb),
                                 "-t", str(cpu)
                                 ]
    print(f"Ejecutando: \n {shlex.join(comando_integronfiltering)}")
    subprocess.run(comando_integronfiltering, 
                                      capture_output = True, text = True, check = True) 
    ruta_genoma = adaptar_dir_integronfiltering(outdir)
    return ruta_genoma


# FUNCIÓN PARA ADAPTAR EL RESULTADO DE INTEGRON FILTERING A INTEGRONFINDER2
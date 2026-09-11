import os
import subprocess 
import shlex 

def adaptar_dir_IntegronFinder(acc_genoma, directorio):
    path_resultados = os.path.join(directorio, f"Results_Integron_Finder_{acc_genoma}")
    resultados_integronfinder = os.listdir(path_resultados)
    for archivo_IntegronFinder in resultados_integronfinder:
        if archivo_IntegronFinder.endswith(".integrons"):
            archivo_pc = os.path.join(path_resultados, archivo_IntegronFinder)
            return archivo_pc
            
def Integron_Finder(genoma, acc_genoma, directorio, args=None): 
    # FIXME -> TENGO QUE REVISAR TODOS LSO ARGS QEU ACEPTA INTEGRONFINDER2
    comando_integronfinder = ["integron_finder", "--promoter-attI",
                                   genoma, "--outdir", directorio]
    print(f"Ejecutando: \n {shlex.join(comando_integronfinder)}")
    subprocess.run(comando_integronfinder, 
                       check=True, capture_output=True, text=True)
    nuevo_archivo = adaptar_dir_IntegronFinder(acc_genoma, directorio)
    return nuevo_archivo


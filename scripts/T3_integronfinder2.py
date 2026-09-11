import os
import subprocess 
import shlex 
from T3_0_integronfiltering import integron_filtering
import pandas as pd

def ajustar_directorio_T2(dir):
    # Procesamos el directorio para que tenga concordancia con el anterior script
    # (antes dabamos -> file = "results/pruebas_T2/kpn.zip")
    archivo = os.path.basename(dir)
    nombre_archivo = os.path.splitext(archivo)[0]
    dir_genomas = os.path.join(os.path.dirname(dir), nombre_archivo, "ncbi_dataset/data")
    return dir_genomas

def acceder_archivos_fasta(dir_in):
    genomas=[]
    carpetas = sorted(os.listdir(dir_in))
    for carpeta in carpetas:
        ruta_individual = os.path.join(dir_in, carpeta)
        if os.path.isdir(ruta_individual):
            archivos_por_carpeta = os.listdir(ruta_individual)
            for archivo in archivos_por_carpeta:
                if os.path.splitext(archivo)[-1] == ".fna":
                    ruta_archivo = ruta_individual + "/" + archivo
                    genomas.append(ruta_archivo)
    return genomas

def adaptar_dir_IntegronFinder(directorio):
    path_resultados=os.path.join(directorio, "Results_Integron_Finder_{}")
    resultados_integronfinder = os.listdir(directorio)
    for carpeta_genoma in resultados_integronfinder:
        path_carpeta_genoma = os.path.join(directorio, carpeta_genoma)
        archivos_genoma = os.listdir(path_carpeta_genoma)
        for archivo in archivos_genoma:
            print(archivo)
            # if archivo.endswith(".integrons"):
            

def Integron_Finder(genoma, directorio, args=None): 
    # FIXME -> TENGO QUE REVISAR TODOS LSO ARGS QEU ACEPTA INTEGRONFINDER2
    comando_integronfinder = ["integron_finder", "--promoter-attI",
                                   genoma, "--outdir", directorio]
    print(f"Ejecutando: \n {shlex.join(comando_integronfinder)}")
    subprocess.run(comando_integronfinder, 
                       check=True, capture_output=True, text=True)
    nuevo_archivo = adaptar_dir_IntegronFinder(directorio)
    

def recorrer_genomas(genomas,
                    dir_out, # FIXME: pongo algo que sea default¿¿¿
                    flag_integronfiltering = True): # Por defecto, se filtra previamente con integron_filtering)
    c=0 ################
    dir_parental = dir_out.replace("/ncbi_dataset/data", "") # como siempre son las msmas carpetas no hay probelma
    for genoma in genomas:
        c+=1 ################
        archivo = os.path.basename(genoma)
        acc = os.path.splitext(archivo)[0]
        if flag_integronfiltering:
            
            directorio_integronfiltering = os.path.join(dir_parental, "results_IntegronFiltering", acc)
            os.makedirs(directorio_integronfiltering, exist_ok=True) 
            # EJECUCIÓN PROGRAMA:
            
            genoma=integron_filtering(genoma, directorio_integronfiltering) # para adaptarlo a la  nueva ruta!
        print(genoma)
        directorio_integronfinder = os.path.join(dir_parental, "results_IntegronFinder2")
        os.makedirs(directorio_integronfinder, exist_ok=True) 
        Integron_Finder(genoma, directorio_integronfinder)
        
        if c==3: break ################

### HAY QUE PONER QUE A INTEGRION SOLO VAA EL ACC PARA HACER EL PAT, PORQUE HACE DOS VECES UN FOR"""!!

    

def extraer_secuencia_Pc(archivo):

    df = pd.read_csv(archivo, 
                     sep="\t",
                     comment="#")

    promotores= df[(df["type_elt"] == "Promoter") & \
                   (df["element"].str.startswith("Pc"))].copy()
    
    print(promotores)

# extraer_secuencia_Pc("/home/marta/CNB/PracticasExternas/results/pruebas_T2/kpn/results_IntegronFinder2/Results_Integron_Finder_GCF_000364385.3_ASM36438v3_genomic.attC_filtered/GCF_000364385.3_ASM36438v3_genomic.attC_filtered.integrons")
    
# EJEMPLOS USADOS:


# path=ajustar_directorio_T2("results/pruebas_T2/kpn.zip")

# genom = acceder_archivos_fasta(path)


# recorrer_genomas(genom, path)
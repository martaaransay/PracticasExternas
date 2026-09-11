import os
import subprocess 
import shlex 
from T3_0_integronfiltering import integron_filtering
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

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



def extraer_secuencia_Pc(ruta_genoma, acc_genoma, archivo_integrons):
    df = pd.read_csv(archivo_integrons, 
                     sep="\t",
                     comment="#")
    #    FIXME -> SE PUEDEN ELIMINAR SIN PORBLEMA OTRAS COLUMNAS Q NO SEAN RELEVANTES...
    promotores= df[(df["type_elt"] == "Promoter") & \
                   (df["element"].str.startswith("Pc"))].copy()
    promotores.insert(0, "ID_genome", acc_genoma)

    genoma_dict = SeqIO.index(ruta_genoma, "fasta")
    
    secuencias = []
    for index, row in promotores.iterrows():
        id_replicon = row['ID_replicon']
        
        # Verificar que el cromosoma/plásmido esté en el archivo FASTA
        if id_replicon not in genoma_dict:
            print(f"Advertencia: {id_replicon} no encontrado en el archivo FASTA. Saltando...")
            continue
            
        # Extraer variables de interés (asegurando formato entero)
        pos_beg = int(row["pos_beg"])
        pos_end = int(row["pos_end"])
        strand = int(row["strand"])
        
        # 3. Extraer la secuencia
        # Python recorta con [inicio:fin]. Restamos 1 a pos_beg para el índice 0.
        secuencia_cruda = genoma_dict[id_replicon].seq[pos_beg - 1 : pos_end]
        
        # 4. Ajustar según la dirección de la hebra (Strand)
        if strand == -1:
            # Obtiene el reverso complementario automáticamente
            secuencia_final = secuencia_cruda.reverse_complement()
        else:
            secuencia_final = secuencia_cruda
            
        # Guardar el resultado
        secuencias.append(str(secuencia_final))
        
    # 5. Añadir las secuencias como una nueva columna al DataFrame original
    df_resultado = promotores.copy()
    df_resultado["sequence"] = secuencias
    
    return df_resultado


def recorrer_genomas(genomas,
                    dir_out, # FIXME: pongo algo que sea default¿¿¿
                    flag_integronfiltering = True): # Por defecto, se filtra previamente con integron_filtering)
    c=0 ################
    dir_parental = dir_out.replace("/ncbi_dataset/data", "") # como siempre son las msmas carpetas no hay probelma
    for genoma in genomas:
        c+=1 ################
        archivo = os.path.basename(genoma)
        acc_genoma = os.path.splitext(os.path.basename(archivo))[0]
        if flag_integronfiltering:         
            directorio_integronfiltering = os.path.join(dir_parental, "results_IntegronFiltering", acc_genoma)
            os.makedirs(directorio_integronfiltering, exist_ok=True) 
            # EJECUCIÓN PROGRAMA:
            genoma=integron_filtering(genoma, directorio_integronfiltering)
            # acutalización del acc (por lo de filtered)
            archivo = os.path.basename(genoma)
            acc_genoma = os.path.splitext(os.path.basename(archivo))[0]
            # para adaptarlo a la  nueva ruta!
        
        directorio_integronfinder = os.path.join(dir_parental, "results_IntegronFinder2")
        os.makedirs(directorio_integronfinder, exist_ok=True) 
        informacion_integronfinder = Integron_Finder(genoma, acc_genoma, directorio_integronfinder)

        ## EXTRAER INFORMACIÓN SOBRE PCs

        pc = extraer_secuencia_Pc(genoma, acc_genoma, informacion_integronfinder)
        print(pc["sequence"])
        if c==3: break ################

### HAY QUE PONER QUE A INTEGRION SOLO VAA EL ACC PARA HACER EL PAT, PORQUE HACE DOS VECES UN FOR"""!!

    


# extraer_secuencia_Pc("/home/marta/CNB/PracticasExternas/results/pruebas_T2/kpn/results_IntegronFinder2/Results_Integron_Finder_GCF_000364385.3_ASM36438v3_genomic.attC_filtered/GCF_000364385.3_ASM36438v3_genomic.attC_filtered.integrons")
    
# EJEMPLOS USADOS:


path=ajustar_directorio_T2("results/pruebas_T2/kpn.zip")

genom = acceder_archivos_fasta(path)


recorrer_genomas(genom, path)
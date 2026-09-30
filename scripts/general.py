import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed 
from pipeline_inicial.T2_downloadgenomes import datasets_flags, download_genomes
from generate_db import download_species_genomes, generate_multifasta, generate_db
from blastn import run_blastn, process_hits, contig_to_genome
from pipeline_inicial.T3_3_Pc import classify_pc_regex, info_to_csv

def create_recovery_file(fetch_directory):
    """docstring"""
    
    fetch_file = os.path.join(fetch_directory, "fetch.txt")
    recovery_file = os.path.join(fetch_directory, "recovery.txt")
    if os.path.isfile(recovery_file):
        pass
    elif os.path.isfile(fetch_file):
        os.replace(fetch_file, recovery_file)  # mueve y sobrescribe
    else:
        print("no fetch file ni recovery file")
        return None
    return recovery_file

def splitting_fetch_file(recovery_file, n_split = 100):
    splits = []
    recovery = []
    with open(recovery_file, "r") as f:
        for line in f:
            splits.append(line)
            if len(splits) == n_split:
                recovery = f.readlines()
                break
    return splits, recovery


def processing(directory, n_split = 100, max_workers = 18, n_tries = 5):
    fetch_directory = os.path.join(directory, "ncbi_dataset")
    recovery_file = create_recovery_file(fetch_directory)
    if recovery_file is None:
        return False
    
    new_fetch_file = os.path.join(fetch_directory, "fetch.txt")
    rehydrate_cmd = ["datasets", "rehydrate", "--directory", directory,
                     "--max-workers", str(max_workers)]

    while True:
        splitted, recovery = splitting_fetch_file(recovery_file, n_split)
        if not splitted:
            break
        with open(new_fetch_file, "w") as f:
            f.writelines(splitted)
        for n in range(1, n_tries+1):
            try:
                subprocess.run(rehydrate_cmd, check=True)
                break
            except subprocess.CalledProcessError as e:
                print(f"Error executing datasets rehydrate:\n{e.stderr}.\nNumber of try: {n}.")       
                if n == n_tries:
                    return False
        with open(recovery_file, "w") as f:
            f.writelines(recovery)
    return True

def process_all(dirs, total_workers=30, n_split=500, n_tries=5):
    # 1) Repartir workers entre especies
    n_especies = len(dirs)                        # 5 especies
    per_dir = max(1, total_workers // n_especies) # 30 // 5 = 6 workers por especie

    failed = []  # aquí apuntaremos las especies que fallen

    # 2) Crear un "equipo" de hilos: uno por especie
    with ThreadPoolExecutor(max_workers=n_especies) as ex:

        # 3) Lanzar una tarea por especie (no espera, solo las encola)
        futs = {}  # diccionario: {tarea_futura: nombre_de_la_especie}
        for name, d in dirs.items():
            fut = ex.submit(processing, d, n_split, per_dir, n_tries)
            futs[fut] = name

        # 4) Ir recogiendo resultados según vayan terminando
        for fut in as_completed(futs):
            name = futs[fut]      # ¿qué especie era?
            ok = fut.result()     # True/False que devolvió processing()

            if ok:
                print(f"{name}: terminado")
            else:
                print(f"{name}: FALLÓ o nada que rehidratar")
                failed.append(name)

    # 5) Al salir del "with", todos los hilos ya han acabado
    return failed

        
if __name__ == "__main__":
    n_genomes = 500
    
    flag_db_created = False
    if not flag_db_created:
        flags = datasets_flags(assembly_level="complete,chromosome,contig",
                            exclude_atypical=True,
                            mag=False)
    
        # species to download
        names = {"Ecoli" : "Escherichia coli", 
                "Klebsiella" : "Klebsiella pneumoniae",
                "Pseudomonas" : "Pseudomonas aeruginosa",
                "Acinetobacter" : "Acinetobacter baumannii",
                "Enterobacter" : "Enterobacter"}

        ## For downloading normal genomes --> download_genomes() directly
            
        file_general = "results/T2_check/EC.zip"
    
        dirs = download_species_genomes(species = names, 
                                    outdir = f"../data/{str(n_genomes)}genomes", 
                                    target_num = n_genomes, 
                                    flags = flags)
        failed = process_all(dirs, total_workers=30, n_split=500)
        for short_name, genome_dir in dirs.items():
            processing(genome_dir)
        if failed:
            print(f"Especies con problemas: {failed}. Relanza el script: recovery.txt retoma donde quedó.")
            exit(1)


    multifasta_path = f"../data/{str(n_genomes)}genomes/multifasta_genomes.fasta"
       
    generate_multifasta(path_to_genomes = f"../data/{str(n_genomes)}genomes/*/*/ncbi_dataset/data/*/*.fna",
                        outfile = multifasta_path)
    generate_db(multifasta = multifasta_path,
                        out_db = f"../data/{str(n_genomes)}genomes/db/db")
    """
    
    query_file = "data/Int1/Int1_Pc_pattern.fa"
    multifasta_path = f"../data/{str(n_genomes)}genomes/multifasta_genomes.fasta"
    genomes_dir = f"../data/{n_genomes}genomes"
    result = run_blastn(query = query_file, 
                        output = f"results/final/{str(n_genomes)}genomes/blast/blast.tsv",
                        db = f"../data/{str(n_genomes)}genomes/db/db",
                        max_target = len(contig_to_genome(genomes_dir))) # max num of seqs to keep: - mapping length

    hits_df = process_hits(result,
                           dir_genomes = f"../data/{n_genomes}genomes",
                           multifasta = multifasta_path)


    for row in hits_df.itertuples():
        pc_info = classify_pc_regex(row.sequence) 
        if pc_info:
            info_to_csv(pc_info, f"results/final/{str(n_genomes)}genomes/final.csv")
"""
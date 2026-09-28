import os
import subprocess
from T2_downloadgenomes import datasets_flags, download_genomes, rehydrate


def create_recovery_file(fetch_directory):
    """docstring"""
    
    fetch_file = os.path.join(fetch_directory, "fetch.txt")
    recovery_file = os.path.join(fetch_directory, "recovery.txt")
    if not os.path.isfile(fetch_file): 
        print("no fetch file")
        return
    else:
        with open(fetch_file, "r") as f_fetch, open(recovery_file, "w") as f_recovery:
            for line in f_fetch:
                f_recovery.write(line)
    rm_cmd = ["rm", fetch_file]
    try:
        subprocess.run(rm_cmd,shell = True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing rm:\n{e.stderr}")
        return 
    return recovery_file


def splitting_fetch_file(recovery_file, n_split=100):
    return


def processing(directory):
    fetch_directory = os.path.join(directory, "ncbi_dataset")
    recovery_file = create_recovery_file(fetch_directory)
    new_fetch_file = os.path.join(fetch_directory, "fetch.txt")
    with open(recovery_file, "r") as file:
        for n_line in range(0,)


if __name__ == "__main__":
    flags_general = datasets_flags()

    file_general = "results/T2_check/EC.zip"

    dir = download_genomes(source = "taxon", 
                     name = "Klebsiella pneumoniae", 
                     filename = file_general, 
                     flags = flags_general,
                     target_num = 250)
    ## FIXME: paralelizar aquí
    splitting_fetch_file
    rehydrate(dir)


recorrer_genomas(file_general)
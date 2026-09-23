import glob
import os
import shlex
import subprocess
import shutil
from pipeline_inicial.T2_downloadgenomes import datasets_flags, download_genomes

def build_database(species, outdir, target_num, flags):
    os.makedirs(outdir, exist_ok = True)
    for short_name, long_name in species.items():
        path = os.path.join(outdir, short_name, f"{short_name}.zip")
        print(path)
        download_genomes(source = "taxon",
                        name = long_name,
                        quality_filter_flag = True,
                        filename = path,
                        flags = flags,
                        target_num = target_num
                        )
        print(f"{short_name} -> download complete")


# Generate multifasta
def generate_multifasta(path_to_genomes, outfile):
    files_to_db = sorted(glob.glob(path_to_genomes))
    with open(outfile, "w") as out:
        for file in files_to_db:
            with open(file, "r") as each_file:
                out.write(each_file.read())
                out.write("\n")

def generate_db(multifasta, out_db):
    # Generate db
    os.makedirs("../data/genomes/db", exist_ok=True)
    makeblastdb_cmd = ["makeblastdb", "-in", multifasta, "-dbtype", "nucl", 
                       "-out", out_db, "-parse_seqids"]
    print(f"Executing: \n{shlex.join(makeblastdb_cmd)}")
    try:
        subprocess.run(makeblastdb_cmd, 
                    check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing makeblastdb:\n{e.stderr}")
        return

def clean_db(outdir):
    #Removes multifasta + genomes + db
    shutil.rmtree(outdir)
    print("Everything removed")
    return

if __name__ == "__main__":
    clean_flag = True
    if clean_flag:
            clean_db(outdir = "data/genomes")

    # # Downloading genomes
    # flags = datasets_flags(assembly_level="complete,chromosome,contig",
    #                        exclude_atypical=True,
    #                        mag=False)
    
    # names = {"Ecoli" : "Escherichia coli", 
    #          "Klebsiella" : "Klebsiella pneumoniae",
    #          "Pseudomonas" : "Pseudomonas aeruginosa",
    #          "Acinetobacter" : "Acinetobacter baumannii",
    #          "Enterobacter" : "Enterobacter spp."}
    
    # build_database(species = names, 
    #                outdir = "data/genomes", 
    #                target_num = 15, 
    #                flags = flags)
    # # Executed once
    # generate_multifasta(path_to_genomes = "data/genomes/*/*/ncbi_dataset/data/*/*.fna",
    #                     outfile = "data/genomes/multifasta_genomes.fasta")
    # generate_db(multifasta = "data/genomes/multifasta_genomes.fasta",
    #             out_db = "data/genomes/db/db")
    
    

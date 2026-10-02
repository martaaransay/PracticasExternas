
from pipeline.downloadgenomes import datasets_flags, download_genomes
from pipeline.generate_db import download_species_genomes, generate_multifasta, generate_db
from pipeline.useful_funct import contig_to_genome
from pipeline.parallel_download import process_all
from pipeline.blastn import run_blastn, process_hits
from pipeline.analyze_Pc import classify_pc_regex, info_to_csv


        
if __name__ == "__main__":
    n_genomes = 10
    
    db_created = False
    if not db_created:
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
            
        dirs = download_species_genomes(species = names, 
                                    outdir = f"../data/{str(n_genomes)}genomes", 
                                    target_num = n_genomes, 
                                    flags = flags)
        failed = process_all(dirs, total_workers=30, n_split=500)
        if failed: 
            print(f"Failed to rehydrate: {failed}")
            exit(0)

    multifasta_path = f"../data/{str(n_genomes)}genomes/multifasta_genomes.fasta"
       
    generate_multifasta(path_to_genomes = f"../data/{str(n_genomes)}genomes/*/*/ncbi_dataset/data/*/*.fna",
                        outfile = multifasta_path)
    generate_db(multifasta = multifasta_path,
                        out_db = f"../data/{str(n_genomes)}genomes/db/db")

    
    query_file = "data/Int1/Int1_Pc_pattern.fa"
    multifasta_path = f"../data/{str(n_genomes)}genomes/multifasta_genomes.fasta"
    genomes_dir = f"../data/{n_genomes}genomes"

    ### FIXME: CAMBIAR LAS COLUMNS PARA AÑADIR DESDE EL JSOn!
    extra_columns = {"bioproject": "assemblyInfo.bioprojectLineage.bioprojects.accession",
                         "biosample": "assemblyInfo.biosample.accession",
                         "organism_name": "organism.organismName",
                         "assembly_level": "assemblyInfo.assemblyLevel"}
    
    result, t = run_blastn(query = query_file, 
                        output = f"results/final/{str(n_genomes)}genomes/blast/blast.tsv",
                        db = f"../data/{str(n_genomes)}genomes/db/db",
                        max_target = len(contig_to_genome(genomes_dir))) # max num of seqs to keep: - mapping length

    hits_df = process_hits(result,
                           dir_genomes = f"../data/{n_genomes}genomes",
                           multifasta = multifasta_path,
                           extra_columns = extra_columns)



    for hit in hits_df.to_dict("records"):
        pc_info = classify_pc_regex(hit["sequence"])

 
        row_info = {"genome": hit["genome"],   
                    "contig": hit["seq_id"],  
                    "species": hit["species"]}
        row_info.update({col: hit[col] for col in extra_columns})
        row_info.update(pc_info)
        info_to_csv(row_info, f"results/final/{str(n_genomes)}genomes/final.csv")



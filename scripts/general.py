from T2_downloadgenomes import datasets_flags, download_genomes, unzip_rehydrate
from T3_complete import
flags_general = datasets_flags()

file_general = "results/T2_check/EC.zip"

check_general = download_genomes(source = "taxon", 
                                 name = "Escherichia coli", 
                                 filename = file_general, 
                                 flags = flags_general)

if check_general:
    unzip_rehydrate(filename = file_general)

recorrer_genomas(file_general)
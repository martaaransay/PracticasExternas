from T2_downloadgenomes import datasets_flags, download_genomes, unzip_rehydrate
from T3_complete import

flags_general = datasets_flags()

file_general = "results/T2_check/EC.zip"

download_genomes(source = "taxon", 
                 name = "Escherichia coli", 
                 filename = file_general, 
                 quality_filter_flag = True,
                 flags = flags_general)

recorrer_genomas(file_general)
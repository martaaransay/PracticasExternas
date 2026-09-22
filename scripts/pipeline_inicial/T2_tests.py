import subprocess
import shlex

from T2_downloadgenomes import datasets_flags, taxon_acc_flag, \
                               quality_filter_genomes, download_genomes, \
                               unzip_rehydrate

#------------------------ Function to define flags for genome download (datasets download)
flags01 = datasets_flags(assembly_source = "refseq",
                        assembly_level = "chromosome,complete",
                        exclude_atypical = True, 
                        exclude_multi_isolate = True
                        )
print(flags01)

flags02 = datasets_flags(assembly_source = "genbank",
                         assembly_level = "complete",
                         mag = "only"
                         )
print(flags02)

flags03 = datasets_flags(annotated = True,
                         exclude_atypical = True, 
                         exclude_multi_isolate = True, 
                         mag = True
                         )
print(flags03)

#------------------------ Function to define taxon/accesion 
# Check not valid arguments -> raise ValueError
# source_name01 = taxon_acc_flag(source = "aa",
#                                name = "Escherichia coli")
# print(source_name01)    

# Check valid arguments
source_name02 = taxon_acc_flag(source = "taxon",
                               name = "Escherichia coli")      
print(source_name02)       

source_name03 = taxon_acc_flag(source = "accession",
                               name = "GCF_000001405.40") 
print(source_name03)    

# Check valid arguments with inputfile -> raise FileNotFoundError
# source_name04 = taxon_acc_flag(source = "inputfile",
#                                name = "prueba.txt") 
# print(source_name04)    


# Check for multiple names (comma-separated)
source_name05 = taxon_acc_flag(source = "taxon",
                               name = "Escherichia coli,Klebsiella pneumoniae") 
print(source_name05)    

# # Control varios elementos en nombre (accesion) modo string
source_name06 = taxon_acc_flag(source = "accession",
                               name = "GCF_000001405.40,GCF_000001635.26") 
print(source_name06)    


##------------------------ Function to filter genomes
##------------------------ Function to download genomes

# Check changing directories

filename01 = "results/T2_check/download01.zip"

# check01 = download_genomes(source = "taxon",
#                            name = "Ginkgo",
#                            quality_filter_flag = True,
#                            filename = filename01,
#                            flags = flags03)

# Check same directory

filename02 = "download02.zip"

# check02 = download_genomes(source = "accession",
#                            name = "GCF_000001405.40,GCF_000001635.26",
#                            quality_filter_flag = True,
#                            filename = filename02)

# Check removing directories
# check03 = download_genomes(source = "taxon",
#                            name = "Ginkgo",
#                            filename = filename02,
#                            flags = flags03)

# # Check removing directories with acc list 
# check04 = download_genomes(source = "taxon",
#                            name = "Ginkgo",
#                            filename = filename02,
#                            quality_filter_flag = True,
#                            flags = flags03)

# # ##------------------------ Function to unzip and rehydrate the genomes
# # if check01:
# #     unzip_rehydrate(filename = filename01)

# # if check02:
# #     unzip_rehydrate(filename = filename02)

# # if check03:
# #     unzip_rehydrate(filename = filename03)

# # if check04:
# #     unzip_rehydrate(filename = filename04)
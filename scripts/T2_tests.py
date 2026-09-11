import subprocess
import shlex 
from T2_descarga_genomas import definir_flags, origen_datos, filtrado_genomas, descargar_genomas, unzip_rehydrate

# ----------------------------------- Función definir_flags()

# Control mayúsculas + False
flags01 = definir_flags(assembly_source="REFSEQ",
                        assembly_level=["chromosome","contIG"],
                        annotated=False,
                        exclude_atypical=False,
                        exclude_multi_isolate=False,
                        mag=False)
# print(flags01)

# Control argumento inválido assemvly_source -> raise ValueError
# flags02 = definir_flags(assembly_source="a",
#                         assembly_level=["chromosome","contIG"],
#                         annotated=False,
#                         exclude_atypical=False,
#                         exclude_multi_isolate=False,
#                         mag=False)

# Control argumento inválido assembly_source -> raise ValueError
# flags03 = definir_flags(assembly_level=["a","contIG"],
#                         annotated=False,
#                         exclude_atypical=False,
#                         exclude_multi_isolate=False,
#                         mag=False)

# Control argumento inválido MAG -> raise ValueError
# flags04 = definir_flags(assembly_source="REFSEQ",
#                         assembly_level=["chromosome","contIG"],
#                         annotated=False,
#                         exclude_atypical=False,
#                         exclude_multi_isolate=False,
#                         mag="O")

# Control assembly_level tipo lista + bool/string/int... para el resto de argumentos
flags05 = definir_flags(assembly_source="REFSEQ",
                        assembly_level="chromosome,contIG, scaffold",
                        annotated=True,
                        exclude_atypical="a",
                        exclude_multi_isolate=1,
                        mag="exclude")
# print(flags05)

# ----------------------------------- Función origen_datos()

# Control origen argumento no válido -> ValueError
# origen_datos01 = origen_datos(origen="aa",
#                               nombre="Escherichia coli")


# Control argumentos válidos con taxon
origen_datos02 = origen_datos(origen="taxon",
                              nombre="Escherichia coli")      
# print(origen_datos02)       

# Control argumentos válidos con accession
origen_datos03 = origen_datos(origen="accession",
                              nombre="GCF_000001405.40") 
# print(origen_datos03)

# Control argumentos válidos con inputfile -> file no encontrado
# origen_datos04 = origen_datos(origen="inputfile",
#                               nombre="prueba.txt") 
# print(origen_datos04)

# Control mayúsculas
origen_datos05 = origen_datos(origen="TAXON",
                              nombre = "Escherichia coli")
# print(origen_datos05)

# Control varios elementos en nombre (taxon) modo string
origen_datos06 = origen_datos(origen="taxon",
                              nombre = "Escherichia coli,Klebsiella pneumoniae")                         
# print(origen_datos06)

# Control varios elementos en nombre (taxon) modo lista
origen_datos07 = origen_datos(origen="taxon",
                              nombre = ["Escherichia coli","Klebsiella pneumoniae"])                         
# print(origen_datos07)

# Control varios elementos en nombre (accesion) modo string
origen_datos08 = origen_datos(origen="accession",
                              nombre = "GCF_000001405.40,GCF_000001635.26")                         
# print(origen_datos08)

# Control varios elementos en nombre (accesion) modo lista
origen_datos09 = origen_datos(origen="accession",
                              nombre = ["GCF_000001405.40","GCF_000001635.26"])                         
# print(origen_datos09)



###-------------------- PARTE FINAL PARA COMPROBAR ERRORES en el comando final

# comando_download = ["datasets", "download", "genome", *filtrado01, 
#                         "--dehydrated", *flags05]

# print(f"Ejecutando: \n {shlex.join(comando_download)}")
# try:
#     subprocess.run(comando_download, check=True, capture_output=True, text=True) 
#     print("TEST PASADO")
# except subprocess.CalledProcessError as e:
#     print(f"Error al ejecutar datasets download:\n{e.stderr}")


# ----------------------------------- Función filtrado_genomas()

# filtrado01 = filtrado_genomas(origen="accession",
#                               nombre = ["GCF_000001405.40","GCF_000001635.26"])
# # print(filtrado01)

# filtrado02 = filtrado_genomas(origen="taxon",
#                               nombre = "Ginkgo")
# # print(filtrado02)


# ----------------------------------- Función descargar_genomas()

# check01=descargar_genomas(origen="accession",
#                   nombre = ["GCF_000001405.40","GCF_000001635.26"],
#                   filename="descarga01.zip",
#                   flags=flags05)


# check02=descargar_genomas(origen="taxon",
#                   nombre = "Ginkgo",
#                   flag_filtrado_calidad=True,
#                   filename="descarga02.zip")


# Control otros directorios

filename03="results/pruebas_T2/descarga03.zip"

# check03 = descargar_genomas(origen="taxon",
#                   nombre = "Ginkgo",
#                   flag_filtrado_calidad=True,
#                   filename=filename03,
#                   flags=flags05)

# Control mismo directorio

# filename04="prueba04.zip"

# check04=descargar_genomas(origen="accession",
#                   nombre = "GCF_000001405.40,GCF_000001635.26",
#                   flag_filtrado_calidad=True,
#                   filename=filename04)

# Control eliminación de carpetas
check05 = descargar_genomas(origen="taxon",
                  nombre = "Ginkgo",
                  filename=filename03,
                  flags=flags05)

# Control eliminación de carpetas con lista de acc (no se elimina)
check06 = descargar_genomas(origen="taxon",
                  nombre = "Ginkgo",
                  filename=filename03,
                  flag_filtrado_calidad=True,
                  flags=flags05)
# ----------------------------------- Función unzip_rehydrate()
# if check03:
#     unzip_rehydrate(filename=filename03)

# if check04:
#     unzip_rehydrate(filename=filename04)

# if check05:
#     unzip_rehydrate(filename=filename03)

if check06:
    unzip_rehydrate(filename=filename03)






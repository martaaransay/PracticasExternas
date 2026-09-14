import subprocess
import shlex 
from T3_0_integronfiltering import integron_filtering

# ----------------------------------- Función integron_filtering()

# EJEMPLO DE GENOMA:
ruta_genoma = "/home/marta/CNB/PracticasExternas/results/pruebas_T2/kpn/ncbi_dataset/data/GCF_000364385.3/GCF_000364385.3_ASM36438v3_genomic.fna"

outdir = "/home/marta/CNB/PracticasExternas/results/pruebas_T2"

a = integron_filtering(ruta_genoma, outdir)

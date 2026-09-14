from Bio import SeqIO

record = SeqIO.read("/home/marta/CNB/PracticasExternas/results/pruebas_T2/kpn/results_IntegronFinder2/Results_Integron_Finder_GCF_000364385.3_ASM36438v3_genomic.attC_filtered/NZ_CP006661.1.gbk", 
                    "genbank")


promotores_pc = []
locus_id = record.id
for feature in record.features:
    if feature.type == "Promoter":
        nombre = feature.qualifiers.get("Promoter", [""])[0]
        if "Pc" in nombre:  # filtra solo Pc, no P_intI1
            secuencia_pc = feature.location.extract(record.seq)
            promotores_pc.append({
                "locus": locus_id,
                "nombre": nombre,
                "coordenadas": str(feature.location),
                "hebra": feature.location.strand,
                "secuencia": str(secuencia_pc),
            })

for p in promotores_pc:
    print(p)
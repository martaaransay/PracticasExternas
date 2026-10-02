# Comparison of BLAST hits vs IF2 integrases in the genomes where they disagree
import os
import re
import glob
import pandas as pd

from pipeline.blastn import process_blast_tsv
from pipeline.useful_funct import contig_to_genome
from pipeline.analyze_Pc import extract_integron_info_gbk

# -------- Paths (copy them from your analyze_methodscomparison.py if they differ)
n = 5                                                    # genomes per species
genomes_dir = f"../data/{n}genomes"
tsv = f"results/{n}genomes/blast/blast_Int1_Pc_pattern.fa.tsv"
if2_dir = f"results/{n}genomes/IF2"
csv = f"results/plots/{n}genomes/genomes_hits.csv"

def genome_key(name):
    match = re.search(r"GC[FA]_(\d+\.\d+)", str(name))
    return match.group(1) if match else str(name)

# -------- 1) Discordant genomes (genome as str to keep the leading zeros)
df = pd.read_csv(csv, dtype = {"genome": str})
disc = df[df["BLAST_hits"] != df["IF2_hits"]]
print("Discordant genomes:")
print(disc[["genome", "BLAST_hits", "IF2_hits"]].to_string(index = False))
disc_genomes = disc["genome"].tolist()

# -------- 2) BLAST rows of those genomes
blast_df = process_blast_tsv(tsv)                        # adds strand, start and end
mapping = contig_to_genome(genomes_dir)
blast_df["genome"] = blast_df["seq_id"].map(
    lambda c: genome_key(mapping[c]) if c in mapping else None)

print(f"\nBLAST rows without genome: {blast_df['genome'].isna().sum()} of {len(blast_df)}")

cols = ["genome", "seq_id", "identity", "length", "query_start", "query_end",
        "start", "end", "strand", "Evalue"]
sub = blast_df[blast_df["genome"].isin(disc_genomes)]
print("\nBLAST rows of the discordant genomes:")
print(sub[cols].sort_values(["genome", "seq_id", "start"]).to_string(index = False))

# -------- 3) IF2 integrases of those genomes (0-based start, as in the gbk)
print("\nIF2 integrases of the discordant genomes:")
for g in disc_genomes:
    gbk_files = []
    for folder in glob.glob(os.path.join(if2_dir, f"Results_Integron_Finder_*{g}*")):
        gbk_files += glob.glob(os.path.join(folder, "*.gbk"))
    info = extract_integron_info_gbk(g, gbk_files)
    for d in info:
        if not d["calin"]:
            print(f"{g}  {d['locus_id']}  {d['start']}-{d['end']}  strand {d['strand']}  "
                  f"{d['integron_type']}  {d['integron_id']}")
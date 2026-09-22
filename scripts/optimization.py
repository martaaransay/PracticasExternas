import subprocess
import shlex
import time
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def run_blastn(query, output, db = "data/genomes/db/db"):
    blast_cmd = ["blastn", "-query", query, "-db", db,
                 "-out", output, "-outfmt", "6"] # añadri evalues¿?¿¿¿
    start = time.time()
    print(f"Executing: \n{shlex.join(blast_cmd)}")
    try:
        subprocess.run(blast_cmd, 
                    check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing blastn:\n{e.stderr}")
        return
    elapsed = time.time() - start
    
    return output

def plot_evalue_distribution(tsv_files):
    columns = ["query_id", "seq_id", "ident", "length", "mismatch", "gapopen",
               "query_start", "query_end", "seq_start", "seq_end", "Evalue", "BitScore"]
    data = []
    log_data = []
    labels = []
    for query, path in tsv_files.items():
        print(query)
        print(path)
        df = pd.read_csv(path, sep="\t", names=columns)
        Evalue = df["Evalue"]
        print(Evalue)
        # data.append(Evalue)
        log_data.append(-np.log10(Evalue)) 
        labels.append(query)
    
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot(log_data, labels=labels)
    ax.set_ylabel("-log10(E-value)")
    ax.set_title("Distribución")
    plt.xticks(rotation=20)
    plt.tight_layout()

    plt.show()
    

# def plot_evalue_distribution(all_hits):
#     fig, ax = plt.subplots(figsize=(7, 5))
#     ax.boxplot(data, labels=labels)
#     ax.set_ylabel("-log10(E-value)")
#     ax.set_title("Distribución de E-values")
#     plt.xticks(rotation=30)
#     plt.tight_layout()

#     plt.show()
#     print(f"[INFO] Plot guardado en {out_path}")


if __name__ == "__main__":
    
    fasta_files = ["data/Int1/Int1.fa", 
                   "data/Int1/Int1_Pc_pattern.fa", 
                   "data/Int1/Int1_no_Pc.fa"]
    all_results = {}
    for each_query in fasta_files:
        name = os.path.basename(each_query)

        each_result = run_blastn(query = each_query, 
                             output = f"results/prueba{name}.tsv")
        all_results[each_query] = each_result
    print(all_results)
    plot_evalue_distribution(all_results)
    # plot_evalue_distribution()

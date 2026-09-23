import subprocess
import shlex
import time
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def run_blastn(query, output, db):
    blast_cmd = ["blastn", "-query", query, "-db", db,
                 "-out", output, "-outfmt", "6"] 
    start = time.time()
    print(f"Executing: \n{shlex.join(blast_cmd)}")
    try:
        subprocess.run(blast_cmd, 
                    check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing blastn:\n{e.stderr}")
        return
    elapsed = time.time() - start

    return output, elapsed

def process_tsv(tsv_files):
    columns = ["query_id", "seq_id", "ident", "length", "mismatch", "gapopen",
               "query_start", "query_end", "seq_start", "seq_end", "Evalue", "BitScore"]
    log_data = []
    bitscore_data = []
    labels = []
    dfs = {}

    for query, path in tsv_files.items():
        df = pd.read_csv(path, sep="\t", names=columns)
        dfs[query] = df
        evalues = df["Evalue"] 
        nonzero = evalues[evalues>0]
        if nonzero.empty:
            min_value = 1e-300
        else:
            min_value = nonzero.min()/10
        evalues = df["Evalue"].replace(0, min_value)  # Replace 0

        log_data.append(-np.log10(evalues)) 
        bitscore_data.append(df["BitScore"])
        labels.append(query)
    return dfs, log_data, bitscore_data, labels

def plot_evalue(logdata, labels, title = None,
                outdir = "results/plots", filename="evalue_distribution.png"):
        # fixme añadir dir out para guardar la imagen

    os.makedirs(outdir, exist_ok=True)


    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot(logdata, tick_labels=labels)
    ax.set_ylabel("-log10(E-value)")
    if title:
        ax.set_title(title)
    else:
        ax.set_title("Distribución")
    plt.xticks(rotation=20)
    plt.tight_layout()

    plt.savefig(os.path.join(outdir, filename))
    plt.show()

def plot_bitscore(bitscore_data, labels, title = None,
                  outdir = "results/plots", filename="bitscore_distribution.png"):
    # fixme añadir dir out para guardar la imagen
    os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.boxplot(bitscore_data, tick_labels=labels)
    ax.set_ylabel("BitScore")
    if title:
        ax.set_title(title)
    else:
        ax.set_title("Distribución")
    plt.xticks(rotation=20)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, filename))
    plt.show()

def plot_evalue_thresholds(dfs, times,
                           title = None,
                           outdir="results/plots", filename="evalue_thresholds.png"):
    os.makedirs(outdir, exist_ok=True)

    range_evalues = np.logspace(1,-300,200)
    fig, axes = plt.subplots(1,2,figsize=(13,5))
    col_counter = 0
    dict_colors = {}
    colors = ["#963FB0", "#41B883", "#1A6D70"]
    for query, df in dfs.items():
        dict_colors[query] = colors[col_counter]
        col_counter += 1
        n_hits = []
        n_unique_hits = []
        for evalue in range_evalues: # For each df, all thresholds are evaluated
            n_hits.append((df["Evalue"] <= evalue).sum())

            n_unique_hits.append(df.loc[df["Evalue"] <= evalue, "seq_id"].nunique())

        axes[0].plot(range_evalues, n_hits, label=query, color=dict_colors[query])
        axes[1].plot(range_evalues, n_unique_hits, label=query, color=dict_colors[query])

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("E-value threshold")
        ax.legend()

    axes[0].set_ylabel("Total Hits")
    axes[1].set_ylabel("Unique Hits (seq_id)")

    plt.tight_layout()
    labels_for_legend = []
    for query, time in times.items():
        labels_for_legend.append(f"{query} - {time:.2f} secs")
    plt.legend(labels = labels_for_legend)

    if title:
        plt.suptitle(title)
    plt.savefig(os.path.join(outdir, filename))
    plt.show()

if __name__ == "__main__":

    fasta_files = ["data/Int1/Int1.fa", 
                   "data/Int1/Int1_Pc_pattern.fa", 
                   "data/Int1/Int1_no_Pc.fa"]
    all_results = {}
    all_times = {}
    for each_query in fasta_files:
        name = os.path.basename(each_query)

        each_result, each_time = run_blastn(query = each_query, 
                                 output = f"results/prueba{name}.tsv",
                                 db = "data/1000genomes/db/db")
        
        all_results[each_query] = each_result
        all_times[each_query] = each_time


    dfs, evalues, bitscore, labels = process_tsv(all_results)
    outdir = "results/plots/1000genomes"
    plot_evalue(evalues, labels, title = "E-Values distribution - db: 5000 genomes", outdir = outdir)
    plot_bitscore(bitscore, labels, title = "BitScore distribution - db: 5000 genomes", outdir = outdir)
    plot_evalue_thresholds(dfs, all_times, title = "E-Values thresholds - db: 5000 genomes", outdir = outdir)

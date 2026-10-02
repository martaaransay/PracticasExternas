# Import libraries 

import os
from xml.parsers.expat import errors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib.patches import Patch

# Plotting
#------------------------ Function to plot e-values distribution
def plot_evalue(logdata, labels, 
                title = "E-value Distribution", outdir = "results/plots", 
                filename = "evalue_distribution.png"):
    """
    Saves a boxplot of the -log10 of e-values for each query.
    
    Args
    ------
    logata (list): A list of lists containing the -log10 of e-values for each query.
    labels (list): A list of the names of the queries.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(7, 5)) # Creates axes for the plot
    # Creates the boxplot with the data and labels defined
    ax.boxplot(logdata, tick_labels=labels) 
    # Adjusting labels / title:
    ax.set_ylabel("-log10(E-value)") # Y label
    ax.set_title(title) # Title of the plot
    plt.xticks(rotation=20) # Rotates the x-axis labels
    plt.tight_layout() # Adjusts the layout to avoid overlapping
    # Saves the plot
    plt.savefig(os.path.join(outdir, filename)) 
    plt.close()
    return

#------------------------ Function to plot bitscore distribution
def plot_bitscore(bitscore_data, labels, 
                  title = "BitScore Distribution", outdir = "results/plots", 
                  filename = "bitscore_distribution.png"):
    """
    Saves a boxplot of the bitscores for each query.
        
    Args
    ------
    bitscore_data (list): A list of lists containing the bitscores for each query.
    labels (list): A list of the names of the queries.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5)) # Creates axes for the plot
    # Creates the boxplot with the data and labels defined
    ax.boxplot(bitscore_data, tick_labels=labels)
    # Adjusting labels / title:
    ax.set_ylabel("BitScore") # Y label
    ax.set_title(title) # Title of the plot
    plt.xticks(rotation=20) # Rotates the x-axis labels
    plt.tight_layout() # Adjusts the layout to avoid overlapping
    # Saves the plot
    plt.savefig(os.path.join(outdir, filename))
    plt.close()
    return

#------------------------ Function to plot e-value thresholds
def plot_evalue_thresholds(dfs, times, mapping=None, key_func=None,
                           evalues=np.logspace(1, -300, 200),
                           colors=["#963FB0", "#41B883", "#1A6D70"],
                           title="E-Values thresholds", outdir="results/plots",
                           filename="evalue_thresholds.png"):
    """
    Saves a plot of the number of hits and unique genomes for each query
    at different e-value thresholds.

    Args
    -----
    dfs (dict): query name -> blastn dataframe.
    times (dict): query name -> execution time.
    mapping (dict, optional): contig -> genome. If given, unique hits are
                              counted as genomes instead of contigs.
    key_func (callable, optional): normalizes the genome name (e.g. genome_key).
    ...
    """
    os.makedirs(outdir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for (query, df), color in zip(dfs.items(), colors):
        df = df.copy()
        # Unit for "unique hits": genome if mapping is given, otherwise contig
        if mapping is not None:
            genome = df["seq_id"].map(mapping)
            if key_func is not None:
                genome = genome.map(lambda g: key_func(g) if pd.notna(g) else g)
            df["unit"] = genome.fillna(df["seq_id"])  # unmapped -> keeps the contig
        else:
            df["unit"] = df["seq_id"]

        n_hits, n_unique = [], []
        for evalue in evalues:
            mask = df["Evalue"] <= evalue
            n_hits.append(mask.sum())
            n_unique.append(df.loc[mask, "unit"].nunique())

        label = f"{query} - {times[query]:.2f} secs"
        axes[0].plot(evalues, n_hits, label=label, color=color)
        axes[1].plot(evalues, n_unique, label=label, color=color)

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("E-value threshold")
        ax.legend()
    axes[0].set_ylabel("Total hits")
    axes[1].set_ylabel("Unique genomes" if mapping is not None else "Unique hits (seq_id)")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)

#------------------------ Function to plot the comparison between methods
def plot_hits_comparison(summary, total_times,
                         title = "Method comparison", outdir = "results/plots",
                         filename = "method_comparison.png"):
    """
    Saves a plot comparing the total hits, unique hits, and execution times across three methods.

    Args
    -----
    summary (dict): A dictionary where keys are the names of the methods and
                    values are tuples with the total hits and unique hits.
    total_times (dict): A dictionary where keys are the names of the methods and
                        values are the execution times.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok=True)

    methods = list(summary.keys()) # List of methods
    n = len(methods) # Number of methods
    # Lists to store results 
    totals = []
    uniques = []
    times = []
    # Colors to plot
    list_colors = ["#35C40E", "#F791F4", "#F1BE3B"]
    colors = {}
    
    for i, meth in enumerate(methods): # Store data from each method
        colors[meth] = list_colors[i]# Set a color for each method
        totals.append(summary[meth][0])
        uniques.append(summary[meth][1])
        times.append(total_times[meth])
        

    fig, axes = plt.subplots(1,2,figsize=(13,5)) # Creates axes for the plot
    x = np.arange(2) # Two comparisons: total and unique
    width = 0.8 / n # Width of each bar

    for i, m in enumerate(methods):
        values = [totals[i], uniques[i]] # Pairs to plot for each method
        center = (n - 1) / 2 # Middle bar 
        offset = width * (i - center) # Horizontal shift
        # Plot the bars
        axes[0].bar(x + offset, # Shift
                    values, width = width, 
                    label = m, color = colors[m])
    # Labels, titles and axes
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(["Total hits", "Unique hits"])
    axes[0].set_ylabel("Number of hits")
    axes[0].set_title("Total hits vs. Unique hits")
    axes[0].legend(title = "Method")

    # Subplot 2: "Time of execution"
    # Background plot to plot times
    axes[1].plot(range(n), times, "-", color = "gray")
    for i, m in enumerate(methods):
        # Plot time of each method
        axes[1].plot(i, times[i], "o", color = colors[m])
        # Annotate the time of each method
        axes[1].annotate(f"{times[i]:.2f}s", (i, times[i]), 
                         textcoords = "offset points", xytext = (0, 10),
                          ha="center")
    
    # Labels, titles and axes
    axes[1].set_xticks(range(n))
    axes[1].set_xticklabels(methods)
    axes[1].set_ylabel("Time (secs)")
    axes[1].set_title("Time of execution")
    axes[1].margins(y=0.15)  # deja hueco arriba para las anotaciones

    plt.suptitle(title) # Set the principal title

    plt.tight_layout()# Adjusts the layout to avoid overlapping
    # Saves the plot
    plt.savefig(os.path.join(outdir, filename))
    plt.close()

    return
    

#------------------------ Function to analyze the overlapping
def plot_hits_overlap(sets_dict, title="Hits overlapping", outdir="results/plots",
                      filename="hits_overlap.png"):
    os.makedirs(outdir, exist_ok=True)
    names = list(sets_dict.keys())
    sets = list(sets_dict.values())
    colors = ["#35C40E", "#F791F4", "#F1BE3B"]   # mismos que plot_hits_comparison

    fig, ax = plt.subplots(figsize=(8, 8))
    venn3(sets, set_labels=("", "", ""), set_colors=colors, alpha=0.5, ax=ax)

    handles = [Patch(facecolor=c, alpha=0.5, label=n) for c, n in zip(colors, names)]
    ax.legend(handles=handles, loc="upper right")

    ax.set_title(title)
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)
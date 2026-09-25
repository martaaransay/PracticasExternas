# Import libraries 
import subprocess
import shlex
import time
import os
from xml.parsers.expat import errors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

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
def plot_evalue_thresholds(dfs, times, evalues = np.logspace(1,-300,200),
                           title = "E-Values thresholds", outdir = "results/plots", 
                           filename = "evalue_thresholds.png"):
    """
    Saves a plot of the number of hits and unique hits for each query at different e-value thresholds.

    Args
    -----
    dfs (dict): A dictionary where keys are the names of the queries and
                values are the corresponding pandas dataframes.
    times (dict): A dictionary where keys are the names of the queries and
                  values are the corresponding execution times.
    evalues (numpy.ndarray, optional): An array of e-values to evaluate.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok=True)
    # Range of e-values to evaluate
    range_evalues = evalues
    # Divides figure into two subplots
    fig, axes = plt.subplots(1,2,figsize=(13,5)) # Creates axes for the plot
    col_counter = 0 # Counter to assign colors to each query
    dict_colors = {} # Dictionary to store colors for each query
    colors = ["#963FB0", "#41B883", "#1A6D70"] # Colors for the queries
    for query, df in dfs.items(): # For each query and its df
        dict_colors[query] = colors[col_counter] # Assigns a color to the corresponding query
        col_counter += 1 # Increase the counter
        # Lists to store hits
        n_hits = [] 
        n_unique_hits = []
        for evalue in range_evalues: # For each df, all thresholds are evaluated
            evalues_in_threshold = df["Evalue"] <= evalue # Stores the evalues <= threshold
            n_hits.append(evalues_in_threshold.sum()) # Counts the hits with <= evalue

            seq_ids_in_threshold = df.loc[evalues_in_threshold, "seq_id"] # Extracts the seq_id column
            n_unique_hits.append(seq_ids_in_threshold.nunique()) # Counts unique hits

        # Plots the data of each evalue threshold
        axes[0].plot(range_evalues, n_hits, label = query, color = dict_colors[query])
        axes[1].plot(range_evalues, n_unique_hits, label = query, color = dict_colors[query])

    # Axis, legends and labels:

    for ax in axes: 
        ax.set_xscale("log") # Log x-axis for the plots
        ax.set_xlabel("E-value threshold") # Same labels for x-axis
    # Set y-axis
    axes[0].set_ylabel("Total Hits")
    axes[1].set_ylabel("Unique Hits (seq_id)")

    # List to store the query + its execution time
    labels_for_legend = [] 
    for query, time in times.items():
        labels_for_legend.append(f"{query} - {time:.2f} secs") # Two decimal numbers
    plt.legend(labels = labels_for_legend) # Set legend

    plt.suptitle(title) # Set the principal title

    plt.tight_layout() # Adjusts the layout to avoid overlapping
    # Saves the plot
    plt.savefig(os.path.join(outdir, filename))
    plt.close()
    return

#------------------------ Function to plot the comparision between methods
def plot_hits_comparison(summary, total_times,
                         title="Method comparision", outdir="results/plots",
                         filename="method_comparison.png"):
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
    n = len(methods)
    # Lists to store results 
    totals = []
    uniques = []
    times = []
    list_colors = ["#35C40E", "#F791F4", "#F1BE3B"]
    colors = {}
    
    for i, meth in enumerate(methods): # Store data from each method
        colors[meth] = list_colors[i]
        totals.append(summary[meth][0])
        uniques.append(summary[meth][1])
        times.append(total_times[meth])
        

    fig, axes = plt.subplots(1,2,figsize=(13,5)) # Creates axes for the plot
    x = np.arange(2) # Two comparisions: total and unique
    width = 0.8 / n

    for i, m in enumerate(methods):
        values = [totals[i], uniques[i]]
        offset = i * width - (width * (n - 1) / 2)
        axes[0].bar(x+offset, values, width=width, label=m, color=colors[m])

    axes[0].set_xticks(x)
    axes[0].set_xticklabels(["Total hits", "Unique hits"])
    axes[0].set_ylabel("Number of hits")
    axes[0].set_title("Total hits vs. Unique hits")
    axes[0].legend(title="Method")

    # Subplot 2: "Time of execution"
    # Plots the times of each method
    axes[1].plot(range(n), times, "-", color="gray", alpha=0.5)
    for i, m in enumerate(methods):
        axes[1].plot(i, times[i], "o", color=colors[m])
        axes[1].annotate(f"{times[i]:.2f}s", (i, times[i]),
                          textcoords="offset points", xytext=(0, 10),
                          ha="center", fontsize=9)

    axes[1].set_xticks(range(n))
    axes[1].set_xticklabels(methods)
    axes[1].set_ylabel("Time (secs)")
    axes[1].set_title("Time of execution")
    axes[1].margins(y=0.15)  # deja hueco arriba para las anotaciones

    for ax in axes:
        ax.tick_params(axis="x", rotation=20)

    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, filename))
    plt.close()
    return
    

#------------------------ Function to analyze the overlapping
def plot_hits_overlap(sets_dict,
                      title="Hits overlapping", outdir="results/plots",
                      filename="hits_overlap.png"):
    """
    Saves a Venn diagram showing the overlap of detected hits between three methods.

    Args
    -----
    sets_dict (dict): A dictionary where keys are the names of the methods and values are sets 
                      with the accession of the genomes with at least one hit for that method.
    title (str, optional): The title of the plot.
    outdir (str, optional): The directory where the plot will be saved.
    filename (str, optional): The name of the file where the plot will be saved.
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok=True)
    names = list(sets_dict.keys()) # List of names of the methods
    sets = list(sets_dict.values()) # List of the sets

    fig, ax = plt.subplots(figsize=(8, 8)) # Creates axes for the plot

    # Executes venn3 function to plot the Venn diagram
    venn3(sets, 
          set_labels = names, # Names of the methods
          ax = ax)
 
    plt.title(title)  # Set the principal title
    # Saves the plot
    plt.savefig(os.path.join(outdir, filename))
    plt.close()
    return

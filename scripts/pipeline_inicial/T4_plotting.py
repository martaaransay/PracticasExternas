# plotting script v0

# Import libraries 
import os
import pandas as pd
import matplotlib.pyplot as plt


path = "/home/mba/Documents/AransayMarta/PracticasExternas/results/final/1500genomes/final.csv"
df = pd.read_csv(path) 
df.columns = df.columns.str.strip() 
count_table = df["Pc_variant"].value_counts(dropna=False)
print(count_table)
# en teoría hay 15 pero beuno
colors = ["#984FB3", "#AAAA00", "#FF5733", 
          "#FF8D1A", "#FFE400", "#FF96E8", "#96FFE5", "#4F80B3"]

plt.figure(figsize=(10, 6))
count_table.plot(kind="bar", color = colors, edgecolor="black")
plt.xlabel("Pc variant")
plt.ylabel("count")
plt.title("Distribution of Pc variants - bbdd: 7500 genomes ESKAPEE")
plt.xticks(rotation=45, ha="right") 
plt.tight_layout()
plt.show()
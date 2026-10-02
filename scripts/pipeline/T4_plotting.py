# Import libraries 
import os
import pandas as pd
import matplotlib.pyplot as plt

path = "/home/mba/Documents/AransayMarta/PracticasExternas/results/final/10genomes/final.csv"
df = pd.read_csv(path) 
df.columns = df.columns.str.strip() 

# Crear una tabla cruzada contando las ocurrencias de Pc_variant por cada species
count_table = pd.crosstab(df["Pc_variant"], df["species"])
print(count_table)

# Colores para las distintas especies (puedes añadir más si ESKAPEE tiene más categorías)
colors = ["#984FB3", "#AAAA00", "#FF5733", "#FF8D1A", "#FFE400", "#FF96E8", "#96FFE5", "#4F80B3"]

plt.figure(figsize=(12, 7))

# Generar el gráfico de barras apiladas
count_table.plot(kind="bar", stacked=True, color=colors[:len(count_table.columns)], edgecolor="black", ax=plt.gca())

plt.xlabel("Pc variant")
plt.ylabel("Count")
plt.title("Distribution of Pc variants by Species - bbdd: 7500 genomes ESKAPEE")
plt.xticks(rotation=45, ha="right") 

# Mover la leyenda fuera del gráfico para que no tape las barras
plt.legend(title="Species", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()

import sys
import pandas as pd
import matplotlib.pyplot as plt


def clasificar_pc_por_especie(cdin, cdout = "abundancia_pc_especie.xlsx"):
    df = pd.read_excel(cdin, 
                       skiprows=3) # las tres primeras líneas se corresponden con el encabezado
    df.columns = df.columns.str.strip() # string y sin espacios en blanco

    tabla_conteo = pd.crosstab(df["Organism"], df["Pc"]) 
    tabla_conteo["Total"] = tabla_conteo.sum(axis=1) #Número total de promotores
    tabla_conteo = tabla_conteo.sort_values("Total", ascending=False) # Ordena del que más al que menos

    # ESKAPEE
    eskapee = ["Escherichia coli", 
               "Klebsiella pneumoniae", 
               "Pseudomonas aeruginosa",
               "Acinetobacter baumannii"] # falta enterobacter que luego se agrupa
    cond_eskapee = tabla_conteo.index.isin(eskapee) | tabla_conteo.index.str.startswith("Enterobacter ")
    tabla_conteo_ESKAPEE = tabla_conteo[cond_eskapee].copy()

    filas_enterobacter = tabla_conteo_ESKAPEE.index.str.startswith("Enterobacter ")
    suma_enterobacter = tabla_conteo_ESKAPEE[filas_enterobacter].sum()
    tabla_conteo_ESKAPEE.loc["Enterobacter spp."] = suma_enterobacter


    with pd.ExcelWriter(cdout) as writer:
        tabla_conteo.to_excel(writer, sheet_name="Conteo")
        tabla_conteo_ESKAPEE.to_excel(writer, sheet_name="Conteo_ESKAPEE")

    especies_plot = eskapee + ["Enterobacter spp."]
    print(especies_plot)
    df_plot = tabla_conteo_ESKAPEE.loc[especies_plot].drop(columns=["Total"]) #escogemos solo las 5 especies genéricas


    print(df_plot)
    # Crear la figura y el gráfico
    plt.figure(figsize=(10, 6))
    df_plot.plot(kind='bar', stacked=True, ax=plt.gca(), edgecolor="black")

    plt.title("Distribución de Promotores en especies ESKAPEE")
    plt.xlabel("Especie", )
    plt.ylabel("Conteo")
    plt.xticks(rotation=45, ha="right")
    #movemos leyenda
    plt.legend(title="Promotor (Pc)", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    plt.show()
    return

a=clasificar_pc_por_especie("/home/marta/CNB/PracticasExternas/data/Integrall/Integrons-Analysis-TJ-V23/Integrons-Analysis-TJ-V23_v2.xlsx")

# Import librerías 
import pandas as pd
import matplotlib.pyplot as plt

# Ruta a la tabla de Integrall
path = "/home/marta/CNB/PracticasExternas/data/Integrall/Integrons-Analysis-TJ-V23/Integrons-Analysis-TJ-V23_v2.xlsx"
encabezado = 3 # Número de líneas en el encabezado a skippear

df = pd.read_excel(path, skiprows = encabezado) # Lectura del excel
df.columns = df.columns.str.strip() # Limpieza de columnas

def abundancia_pc(df, 
                  columnas, 
                  path_out = "./results/pruebas_inicial_integrall/abundancia_pc_especie.xlsx"):
    """
    Genera un xlsx de abundancia de promotores por especie y un conjunto de gráficos.

    Args
    ----------
    df (pd.DataFrame): Tabla de entrada con dos columnas (especies y Pcs)
    columnas (list): Lista con el nombre de las columnas que contienen las especies y los Pcs.
    path_out (str, optional): Ruta del archivo .xlsx de salida donde se guardarán las tablas de conteo.
    """
    tabla_conteo = pd.crosstab(df[columnas[0]], df[columnas[1]]) # Conteo del número de Pc - especie
    # Columna con el conteo total de Pc por especie
    tabla_conteo["Total"] = tabla_conteo.sum(axis=1) 
    # Ordenación de la tabla de mayor a menor número de Pc
    tabla_conteo = tabla_conteo.sort_values("Total", ascending=False) # Ordena del que más al que menos

    # ------------Filtro especies ESKAPEE------------
    eskapee = ["Escherichia coli", 
                "Klebsiella pneumoniae", 
                "Pseudomonas aeruginosa",
                "Acinetobacter baumannii",
                "Enterobacter spp."] 
    # Condición especie ESKAPEE 
    especies = tabla_conteo.index # Nombre de las filas
    # Condición: Coincida con ESKAPEE o empiece por Enterobacter (.str para aplicar startswith)
    cond_eskapee = especies.isin(eskapee) | especies.str.startswith("Enterobacter ") 

    # Nueva tabla con el conteo único para el filtro
    tabla_conteo_ESKAPEE = tabla_conteo[cond_eskapee].copy()
    # Nueva fila con el sumatorio de todas las Enterobacter
    cond_enterobacter = tabla_conteo_ESKAPEE.index.str.startswith("Enterobacter ")
    suma_enterobacter = tabla_conteo_ESKAPEE[cond_enterobacter].sum()
    # Adición de la nueva fila generada
    tabla_conteo_ESKAPEE.loc["Enterobacter spp."] = suma_enterobacter
 
    # ----------Nuevo xlsx con ambos conteos-----------
    with pd.ExcelWriter(path_out) as w:
        tabla_conteo.to_excel(w, sheet_name="Conteo")
        tabla_conteo_ESKAPEE.to_excel(w, sheet_name="Conteo_ESKAPEE")

    # -----------Gráficas de ambos conteos------------
    # Colores para los Pc:
    colores = ["#A7A8A7", "#B3924F", "#984FB3", "#AAAA00", "#FF5733", 
               "#FF8D1A", "#FFE400", "#FF96E8", "#96FFE5", "#4F80B3" ]
    
    plt.figure(figsize=(10, 6))

    # Gráfica 1: Abundancia general (sin tener en cuenta la especie)
    plt.subplot(1, 2, 1)
    # Eliminación columna Total (no es Pc a printear)
    df_plot = tabla_conteo.drop(columns=["Total"]).T # Transposición de la tabla
    # Nuevo df con la abundancia general de los promotores (sin separar por especie individual)
    df_total = df_plot.sum(axis=1) 

    # Generación del gráfico de barras con los colores definidos anteriormente 
    df_total.plot(kind="bar", color=colores, edgecolor="black")

    # Etiquetas del gráfico
    plt.title("Abundancia relativa de promotores", fontweight="bold")
    plt.xlabel("Promotores", fontweight="bold")
    plt.ylabel("Conteo", fontweight="bold")
    plt.xticks(rotation=45, ha="right")

    # Gráfica 2: Especies ESKAPEE
    plt.subplot(1, 2, 2)
    
    # Filtro de las especies ESKAPEE y eliminación columna Total 
    df_plot_eskapee = tabla_conteo_ESKAPEE.loc[eskapee].drop(columns=["Total"]) 

    df_plot_eskapee.plot(kind="bar", 
                         ax=plt.gca(), # Para que no genere dos gráficos
                         stacked=True, # Barras apiladas
                         color=colores, 
                         edgecolor="black")

    plt.title("Distribución de Promotores en especies ESKAPEE", fontweight="bold")
    plt.xlabel("Especies ESKAPEE", fontweight="bold")
    plt.ylabel("Conteo", fontweight="bold")
    plt.xticks(rotation=45, ha="right")
    # Leyenda de los gráficos (colores de los Pc)
    plt.legend(title="Promotor (Pc)", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout() # AJuste de tamaño
    plt.show()
    return

abundancia_pc(df, ["Organism", "Pc"])

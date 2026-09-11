path = "/home/marta/CNB/PracticasExternas/results/pruebas_T2/kpn/results_IntegronFiltering/GCF_000364385.3_ASM36438v3_genomic/GCF_000364385.3_ASM36438v3_genomic.attC_filtered.fasta"
diccionario_fasta = {} # Es mejor no usar 'dict' como nombre de variable para no sobreescribir la función nativa

with open(path, "r") as archivo:
    id_actual = ""
    secuencia_actual = []
    
    for linea in archivo:
        linea = linea.strip() # Elimina el salto de línea al final
        
        if linea.startswith(">"):
            # Si ya teníamos un ID previo, unimos su secuencia y la guardamos en el diccionario
            if id_actual:
                diccionario_fasta[id_actual] = "".join(secuencia_actual)
            
            # Guardamos el nuevo ID y reiniciamos la lista para su secuencia
            id_actual = linea[1:] 
            secuencia_actual = [] 
        else:
            secuencia_actual.append(linea)
            
    # Procesar el último registro del archivo al salir del bucle
    if id_actual:
        diccionario_fasta[id_actual] = "".join(secuencia_actual)

# Comprobación de que funciona
seq = diccionario_fasta["NZ_CP006661.1"]
print(seq[114955:114981])

seq2 = diccionario_fasta["NZ_CP006662.2"]
print(seq2[8128:8154])
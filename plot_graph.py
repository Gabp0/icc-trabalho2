import os
import string
import matplotlib.pyplot as plt
import pandas as pd

if not os.path.isdir("graficos"):
    os.makedirs("graficos")

sizes = [10, 32, 50, 64, 100, 128, 200, 250, 256, 300, 400, 512, 600, 1000, 1024, 2000, 2048, 3000, 4096]
preOtm_df = pd.read_csv('likwid_output/filtered-pre-otm.csv').T
afterOtm_df = pd.read_csv('likwid_output/filtered-after-otm.csv').T

def plotGraph(metric:list, logy:bool, title:string, metric_name:string):
# gera o grafico e exporta
    metric_df = pd.DataFrame(data=metric)

    ax:plt.Axes = metric_df.plot(marker='.', logy=logy, grid=True, title=title, 
                                 xlabel="Tamanho", ylabel=metric_name) 
    ax.set_xticks(range(len(sizes)))
    ax.set_xticklabels(sizes, rotation=90)

    plt.savefig(f"graficos/{title.replace(' ', '')}.png")
    plt.close()

# graficos de tempo
for i in range(0,7):    # offset das linhas com os valores
    time = []

    title = f"Tempo - Método de{preOtm_df.iloc[1, (i*5) + 1].split('_')[0].replace('Region', ' ')}"
    for j in range(0, len(sizes)): 
        time.insert(j, {"Nao-otimizado" : float(preOtm_df.iloc[1, (40*j + (i*5))]), 
                           "Otimizado" : float(afterOtm_df.iloc[1, (40*j + (i*5))])})
    
    plotGraph(time, True, title, "Templo (Escala Log)")

    title = f"Operações Aritiméticas - Método de{preOtm_df.iloc[1, (i*5) + 1].split('_')[0].replace('Region', ' ')}"
    art_op = []
    for j in range(0, len(sizes)): 
        art_op.insert(j, {"Nao-otimizado DP" : float(preOtm_df.iloc[1, (40*j + (i*5) + 2)]), 
                          "Nao-otimizado AVX" : float(preOtm_df.iloc[1, (40*j + (i*5) + 3)]),
                          "Otimizado DP" : float(afterOtm_df.iloc[1, (40*j + (i*5) + 2)]), 
                          "Otimizado AVX" : float(afterOtm_df.iloc[1, (40*j + (i*5) + 3)])})
        
    plotGraph(art_op, False, title, "Operações Aritiméticas")

    mem_bd = []
    title = f"Banda de memória - Método de{preOtm_df.iloc[1, (i*3) + 1215].split('_')[0].replace('Region', ' ')}"
    for j in range(0, len(sizes)): 
        mem_bd.insert(j, {"Nao-otimizado" : float(preOtm_df.iloc[1, (24*j + (i*3) + 1217)]), 
                           "Otimizado" : float(afterOtm_df.iloc[1, (24*j + (i*3) + 1217)])})

    plotGraph(mem_bd, False, title, "Largura de banda de memória")

    cache_miss = []
    title = f"Taxa de erro de cache - Método de{preOtm_df.iloc[1, (i*3) + 759].split('_')[0].replace('Region', ' ')}"
    for j in range(0, len(sizes)): 
        cache_miss.insert(j, {"Nao-otimizado" : float(preOtm_df.iloc[1, (24*j + (i*3) + 761)]), 
                           "Otimizado" : float(afterOtm_df.iloc[1, (24*j + (i*3) + 761)])})

    plotGraph(cache_miss, False, title, "Taxa de erro de cache")

    
# Gabriel de Oliveira Pontarolo GRR20203895
# Rodrigo Saviam Soffner GRR20205092
# Script para geracao dos gráficos
import os
import string
import matplotlib.pyplot as plt
import pandas as pd

if not os.path.isdir("graficos"):
    os.makedirs("graficos")

sizes = [10, 32, 50, 64, 100, 128, 200, 250, 256, 300, 400, 512, 600, 1000, 1024, 2000, 2048, 3000, 4096] # tamanho das funcoes
preOtm_df = pd.read_csv('likwid_output/filtered-pre-otm.csv').T # nao otimizado
afterOtm_df = pd.read_csv('likwid_output/filtered-after-otm.csv').T # otimizado

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
for i in range(0,8):    # offset das linhas com os valores
    
    # tempo de execucao
    time = []
    title = f"Tempo - Método de{preOtm_df.iloc[1, (i*5) + 1].split('_')[0].replace('Region', ' ')}"
    print(title)
    for j in range(0, len(sizes)):  # separa os valores em uma lista de dicionarios
        time.insert(j, {"Nao-otimizado" : float(preOtm_df.iloc[1, (40*j + (i*5))]), 
                           "Otimizado" : float(afterOtm_df.iloc[1, (40*j + (i*5))])})
    
    plotGraph(time, True, title, "Templo (Escala Log)")

    
    # dp e avx
    title = f"Operações Aritiméticas - Método de{preOtm_df.iloc[1, (i*5) + 1].split('_')[0].replace('Region', ' ')}"
    print(title)
    art_op = []
    for j in range(0, len(sizes)): # separa os valores em uma lista de dicionarios
        art_op.insert(j, {"Nao-otimizado DP" : float(preOtm_df.iloc[1, (40*j + (i*5) + 2)]), 
                          "Nao-otimizado AVX" : float(preOtm_df.iloc[1, (40*j + (i*5) + 3)]),
                          "Otimizado DP" : float(afterOtm_df.iloc[1, (40*j + (i*5) + 2)]), 
                          "Otimizado AVX" : float(afterOtm_df.iloc[1, (40*j + (i*5) + 3)])})
        
    plotGraph(art_op, False, title, "MFlops/s")

    # largura de banda de memoria
    mem_bd = []
    title = f"Banda de memória - Método de{preOtm_df.iloc[1, (i*3) + len(sizes)*64 - 1].split('_')[0].replace('Region', ' ')}"
    print(title)
    for j in range(0, len(sizes)): # separa os valores em uma lista de dicionarios
        mem_bd.insert(j, {"Nao-otimizado" : float(preOtm_df.iloc[1, (24*j + (i*3) + len(sizes)*64 + 1)]), 
                           "Otimizado" : float(afterOtm_df.iloc[1, (24*j + (i*3) + len(sizes)*64 + 1)])})

    plotGraph(mem_bd, False, title, "MFlops/s")

    # taxa de erro da cache
    cache_miss = []
    title = f"Taxa de erro de cache - Método de{preOtm_df.iloc[1, (i*3) + len(sizes)*40 - 1].split('_')[0].replace('Region', ' ')}"
    print(title)
    for j in range(0, len(sizes)): # separa os valores em uma lista de dicionarios
        cache_miss.insert(j, {"Nao-otimizado" : float(preOtm_df.iloc[1, (24*j + (i*3) + len(sizes)*40 + 1)]), 
                           "Otimizado" : float(afterOtm_df.iloc[1, (24*j + (i*3) + len(sizes)*40 + 1)])})

    plotGraph(cache_miss, False, title, "Taxa de erro de cache")

    

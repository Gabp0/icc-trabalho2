#!/bin/bash
# Gabriel de Oliveira Pontarolo GRR20203895
# Rodrigo Saviam Soffner GRR20205092
# Script para execucao dos testes no programa
#uso: ./run.sh {caminho para arquivo de entrada}

METRICA="FLOPS_DP L3 L2CACHE" #L3 L2CACHE ENERGY
CORE=3

# testa o arquivo de entrada
if [ -z "$1" ]
then
    INPUT="rosenbrock.txt"
else
	INPUT=$i
fi

# testa os diretorios de saida
if [ ! -d likwid_output ]
then
    mkdir likwid_output
    echo "Diretório likwid_output/ criado"
fi

if [ ! -d newtonpc_output ]
then
    mkdir newtonpc_output
    echo "Diretorio newtonpc_output/ criado"
fi

echo "Compilando"
make -C pre-otimizacao
make -C otimizado

for m in ${METRICA}
do
    echo "Executando o programa nao otimizado para a metrica $m"
    likwid-perfctr -o likwid_output/pre-otm_${m}.csv  -O -C ${CORE} -g ${m} -m ./pre-otimizacao/newtonPC < $INPUT > newtonpc_output/pre-otm_${m}.txt
    echo "Executando o programa otimizado para a metrica $m"
    likwid-perfctr -o likwid_output/after-otm_${m}.csv  -O -C ${CORE} -g ${m} -m ./otimizado/newtonPC < $INPUT > newtonpc_output/after-otm_${m}.txt
done

# filtro inicial

OTIMIZACOES="pre-otm after-otm"

for otm in ${OTIMIZACOES}
do
    egrep "Group 1 Raw|RDTSC Runtime|DP|AVX DP" likwid_output/${otm}_FLOPS_DP.csv >> likwid_output/filtered-${otm}.csv
    egrep "Group 1 Raw|RDTSC Runtime|L2 miss ratio" likwid_output/${otm}_L2CACHE.csv >> likwid_output/filtered-${otm}.csv
    egrep "Group 1 Raw|RDTSC Runtime|L3 bandwidth" likwid_output/${otm}_L3.csv >> likwid_output/filtered-${otm}.csv
done

# gera os graficos

echo "Gerando os gráficos"

python3 plot_graph.py

echo "Pronto"


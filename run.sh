#!/bin/bash
#uso: ./run.sh {caminho para arquivo de entrada}

METRICA="FLOPS_DP L3 L2CACHE ENERGY" 
CORE=3

if [ -z "$1" ]
then
	echo "Especifique o caminho para a entrada: ./run.sh {caminho para arquivo de entrada}"
	exit
fi


if [ ! -d likwid_output ]
then
	mkdir likwid_output
fi

if [ ! -d newtonpc_output ]
then
	        mkdir newtonpc_output
fi

make -C pre-otimizacao
make -C otimizado

for m in ${METRICA}
do
	likwid-perfctr -o likwid_output/pre-otm_${m}.csv  -O -C ${CORE} -g ${m} -m ./pre-otimizacao/newtonPC < $1 > newtonpc_output/pre-otm_${m}.txt
	likwid-perfctr -o likwid_output/after-otm_${m}.csv  -O -C ${CORE} -g ${m} -m ./otimizado/newtonPC < $1 > newtonpc_output/after-otm_${m}.txt
done

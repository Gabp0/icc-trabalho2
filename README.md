# **Introdução a Computação Científica - Trabalho 2**

## Autores:

- Gabriel de Oliveira Pontarolo GRR20203895
- Rodrigo Saviam Soffner GRR20205092

## Arquivos:

- **run.sh** - Script em bash que compila, executa com o LIKWID, filtra a saída e chama o script em python para geração dos gráficos
- **plot_graph.py** - Script em python que gera os gráficos a partir da saída filtrada do LIKWID
- **rosenbrock.txt** - Arquivo contendo as Funções de Rosenbrock de dimensões entre 10 e 4096
- **RELATORIO-gop20-rss20.pdf** - Arquivo .pdf contendo o relatório sobre os experimentos realizados
- Diretório **otimizado/** - Código fonte da versão otimizada do código
- Diretório **pre-otimizacao/** - Código fonte da versão original do código

## Compilação, execução e geração dos gŕaficos:

    Dentro do diretório **gop20-rss20/**, basta executar *./run.sh [entrada (opcional)]* e o script irá realizar as três etapas sozinho.
    Inicialmente, ele irá compilar os dois códigos fonte utilizando seus respectivos Makefiles, ambos com as flags de otimização *-O3 -mavx -march=native -fstrict-aliasing* e com as flags para o uso do LIKWID.
    Após isso, ele irá executar os dois programas com o *likwid-perfctr*, com seus respectivos parâmetros para cada uma das métricas: **FLOPS_DP**, **L3** e **L2CACHE**, utilizando como entrada o parâmetro opcional passado durante a chamada do script ou o arquivo **rosenbrock.txt** caso nenhuma entrada tenha sido especificada.
    As saídas do programa são direcionadas para o diretório, criado pelo próprio script, **newtonpc_output/** com o formato **pre-otm_[métrica].txt** para a versão não otimizada e **after-otm_[métrica].txt** para a otimizada.
    As saídas do *likwid-perfctr* são direcionadas para o diretório, também criado pelo script, **likwid_output/** com o mesmo formato **pre-otm_[métrica].csv** para a não otimizada e **after-otm_[métrica].csv** para a otimizada.
    Quando a execução dos 6 experimentos termina, é feita uma filtragem inicial, fazendo uso do comando *egrep*, separando os valores de interesse de cada uma das saídas em um único arquivo (*filtered-pre-otm.csv* para o não otimizado e *filtered-after-otm.csv* para o otimizado). Mais especificamente, os valores **RDTSC Runtime**, **DP [MFLOPS/s], **AVX DP [MFLOPS/s]**, **L2 miss ratio** e **L3 bandwidth**.
    Finalmente, é executado o script em python **plot_graph.py** que irá gerar cada um dos 32 gŕaficos no diretório **graficos/** com o formato **[métrica]-[método]-[região].png**, utilizando como base os dois arquivos com a saída filtrada do *likwid-perfctr*.
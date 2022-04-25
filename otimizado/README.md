# **Introdução a Computação Científica - Trabalho 1**

## Autores:

- Gabriel de Oliveira Pontarolo GRR20203895
- Rodrigo Saviam Soffner GRR20205092

## Uso:

- ./newtonPC [-o "saída"]

## Módulos:

  - [**Main.c**](#mainc)
  - [**Functions.h**](#functionsh)
  - [**GaussianElimination.h**](#gaussianeliminationh)
  - [**GaussSeidel.h**](#gaussseidelh)
  - [**LuDecomposition.h**](#ludecompositionh)
  - [**Utils.h**](#utilsh)
    - [**Códigos de Saída de Erro**](#saidas-de-erro)
  - [**Métodos de Newton**](#métodos-de-newton)
    - [**Newton Padrão**](#newton-padrao)
    - [**Newton Inexato**](#newton-inexato)
    - [**Newton Modificado**](#newton-modificado)


### Main.c
  O programa principal consiste de um loop no qual a expressão matemática é lida da entrada padrão, os três métodos de newton são executados e o resultado é impresso no arquivo de saída especificado ou na saída padrão, caso nenhum arquivo tenha sido passado como argumento. Após isso, ele libera a memória utilizada e lê a próxima enquanto a entrada for diferente de 'EOF'.

### Métodos de Newton
  As três variações do método de newton utilizados para calcular o sistema linear gerado pela Hessiana e Gradiente da respectiva expressão matemática. O resultado dos três é armazenado em seu respectivo ponteiro para a struct N_RESULT pela biblioteca ["functions.h"](#functionsh) o que inclui o número de iterações necessárias pelo método, o valor de f de cada iteração e os tempos gastos.
#### Newton Padrão
  Implementado na biblioteca "newtonPadrao.h", utiliza o método da [eliminação de Gauss](#guassianeliminationh) para a resolução do sistema linear.
#### Newton Modificado
  Implementado na biblioteca "newtonModificado.h", utiliza o método da [Fatoração LU](#ludecompositionh) para a resolução do sistema linear.
#### Newton Inexato
  Implementado na biblioteca "newtonInexato.h", utiliza o método iterativo de [Gauss-Seidel](#gaussseidelh) para a resolução do sistema linear.

### Functions.h
  Funções que gerenciam a entrada e saída das expressões matemáticas. Inclue a leitura e a alocação de memória através da função "*FUNCTION *readFunction(void);*" que lê da entrada padrão no formato:
  "Número de variaveis
   Expressão da função
   Aproximação inicial das variáveis
   Número máximo de iterações"
assim como a escrita dos resultados dos métodos com a função *"void printMethod(FUNCTION \*func, char \*output);"* e o cálculo de ambos o Gradiente (*void Gradiente(*FUNCTION *func, void \**grad);*) e a Hessiana (*void Hessiana(FUNCTION *func, void \**grad, void \*\**hessi);*) da expressão. As operações nas expressões são realizadas através da biblioteca [matheval](https://www.gnu.org/software/libmatheval/).
  
### GaussianElimination.h
  Funções para representar e resolver sistemas lineares utilizando o método da eliminação de Gauss. Inclui a alocação de memória com *LINEAR_SYST *initLS(int size);* e a resolução propriamente dita do sistema (*void gaussianElimination(LINEAR_SYST *syst);*).
  
### LuDecomposition.h
  Funções para representar e resolver sistemas lineares utilizando o método da fatoração LU. Ao contrário da [eliminação de Gauss](#gaussianeliminationh), essa biblioteca representa o sistema utilizando as matrizes triangular superior U e triangular inferior L. Inclui a alocação de memória com *LINEAR_SYST_LU *initLSLU(int size);*, a fatoração da matriz de coeficientes em L e U com *void factorize(LINEAR_SYST_LU *syst);* (note que a matriz de coeficientes deve ser passada como argumento para função dentro da matriz de double L** da struct LINEAR_SYST_LU, esta será substituída pela respectiva triangular superior cálculada) e a resolução própriamente dita do sistema (*void solveLU(LINEAR_SYST_LU *syst);*). As matrizes L e U são recalculadas apenas a cada HESS_STEPS iterações, sendo esse o número de variáveis da expressão matemática.
  
### GaussSeidel.h 
 Funções para representar e resolver sistemas lineares utilizando o método iterativo de Gauss-Seidel. A struct utilizada para representar o sistema é similar a da [eliminação de Gauss](#gaussianeliminationh), com o adicional do vetor de double X cálculado na iteração anterior do método, o qual é usado para testar um dos critérios de parada do método. A função *LINEAR_SYST_GS *initLSGS(int size);* faz a alocação de memória e a função *void gaussSeidel(LINEAR_SYST_GS *syst);* resolve o sistema. O número máximo de iterações é 50 e a tolerância é 10^-6.
 
### Utils.h
  Funções simples para uso geral em diversas partes do programa. Em especial *char *getArgs(int argc, char \**argv);* que separa os argumentos lidos na [chamada do programa](#uso) e retorna o nome do arquivo de saída ou NULL caso nenhum argumento tenha sido passado (o programa imprime na saída padrão nesse caso). Há também a função *void exitStatus(int code);* que gerencia a saída de erro.

#### Códigos de saídas de erro
  Códigos de saída de erro do programa gerados pela função *void exitStatus(int code);* da biblioteca [utils.h](#utilsh)
- 101 : Falha na alocação de memória 
- 102 : Ponteiro aponta para uma localidade inválida de memória
- 103 : Caso houve erro em uma das funções da biblioteca [matheval](https://www.gnu.org/software/libmatheval/)
- 104 : Caso o número de argumentos passados pela linha de comando é diferente de dois (caso o arquivo de saída tenha sido especificado) ou nenhum
- 105 : Argumento passado pela linha de comando é diferente do [especificado](#uso).
- 106 : Caso de falha ao abrir um arquivo pela função "fopen()"
- 107 : Caso houve divisão por zero em uma das operações
  Caso não haja nenhum erro, o programa termina sua execução com o código de saída 0;



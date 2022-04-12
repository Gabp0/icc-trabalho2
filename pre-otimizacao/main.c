// Introdução a Computação Científica - Trabalho 1
// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Programa principal

#include "functions.h"
#include "newtonPadrao.h"
#include "newtonModificado.h"
#include "newtonInexato.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    char *output = getArgs(argc, argv); // argumentos pela linha de comando
    FUNCTION *input_func;

    do
    {
        input_func = readFunction(); // le a funcao

        // metodos de resolucao de sistema linear
        NewtonPadrao(input_func);
        NewtonModificado(input_func);
        NewtonInexato(input_func);

        printMethod(input_func, output);

        deleteFunction(input_func); // libera memoria
    } while (getc(stdin) != EOF);
    return EXIT_SUCCESS;
}
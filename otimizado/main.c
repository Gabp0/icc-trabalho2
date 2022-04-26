// Introdução a Computação Científica - Trabalho 2 - Otimizado
// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Programa principal

#include "functions.h"
#include "newtonPadrao.h"
#include "newtonInexato.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <likwid.h>

int main(int argc, char **argv)
{

    char *output = getArgs(argc, argv); // argumentos pela linha de comando
    FUNCTION *input_func;

    LIKWID_MARKER_INIT;

    do
    {
        input_func = readFunction(); // le a funcao

        string_t markerNewtonPadrao = markerName("NewtonPadrao", input_func->var_num);
        string_t markerNewtonInexato = markerName("NewtonInexato", input_func->var_num);

        // metodos de encontrar as raizes da funcao
        LIKWID_MARKER_START(markerNewtonPadrao);
        NewtonPadrao(input_func);
        LIKWID_MARKER_STOP(markerNewtonPadrao);

        LIKWID_MARKER_START(markerNewtonInexato);
        NewtonInexato(input_func);
        LIKWID_MARKER_STOP(markerNewtonInexato);

        printMethod(input_func, output);

        deleteFunction(input_func); // libera memoria
    } while (getc(stdin) != EOF);

    LIKWID_MARKER_CLOSE;

    return EXIT_SUCCESS;
}
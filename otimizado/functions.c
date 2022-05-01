// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacoes das funcoes da biblioteca functions.h

#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

N_RESULT *_initNR(void)
// aloca memoria para a struct de resultados do metodo de newton
{
    N_RESULT *new = malloc(sizeof(N_RESULT));
    if (!new)
        exitStatus(MEM_ALOC);

    new->it_num = 0;

    return new;
}

void _deleteResult(N_RESULT *result)
// libera memoria da struct de resultados do metodo de newton
{
    free(result->f_k);
    free(result);
}

FUNCTION *readFunction(void)
// Le a funcao da entrada padrao e armazena em uma struct do tipo FUNCTION. retorna o ponteiro para a struct
{
    FUNCTION *function = malloc(sizeof(FUNCTION));

    if (!function)
        exitStatus(MEM_ALOC);

    function->var_num = 0;

    fscanf(stdin, "%d\n", &function->var_num);
    fscanf(stdin, "%*[^\n]\n");
    if (function->var_num == 0)
        exit(EXIT_SUCCESS);

    function->initial_aps = malloc(sizeof(double) * function->var_num);
    // function->names = malloc(sizeof(char **) * function->var_num);
    if (!function->initial_aps) // || (!function->names))
        exitStatus(MEM_ALOC);

    for (int i = 0; i < function->var_num; i++) // armazena os nomes das variaveis
        fscanf(stdin, "%lf", &function->initial_aps[i]);

    fscanf(stdin, "%lf\n%d", &function->t_ep, &function->it_num); // epsilon e numero maximo de iteracoes

    // aloca memoria para armazenar os resultados dos 3 metodos de newton
    function->n_p = _initNR();
    function->n_m = _initNR();
    function->n_i = _initNR();

    return function;
}

void printMethod(FUNCTION *func, char *output)
// imprime o resultado dos 3 metodos de newton e seus tempos de execucao para a funcao func
{
    FILE *output_file; // arquivo de saida onde sera impresso
    if (!output)
        output_file = stdout;
    else
        output_file = fopen(output, "a");

    if (!output_file)
        exitStatus(FOPEN_ERR);

    // cabeçalho
    fprintf(output_file, "%d\n#", func->var_num); //, func->expression);
    fprintf(output_file, "Iteração \t| Newton Padrão \t| Newton Inexato\n");

    int z = max(func->n_p->it_num, func->n_m->it_num, func->n_i->it_num); // encontra o metodo com o maior numero de iteracoes
    for (int i = 0; i < z; i++)
    {
        // imprime iteração
        fprintf(output_file, "%d \t\t| ", i);

        if (func->n_p->it_num > i)
            fprintf(output_file, "%1.14e\t| ", func->n_p->f_k[i]); // newton padrao
        else
            fprintf(output_file, "\t\t\t| ");

        if (func->n_i->it_num > i)
            fprintf(output_file, "%1.14e\t \n", func->n_i->f_k[i]); // newton inexato
        else
            fprintf(output_file, "\t\t\t \n");
    }

    if (output)
        fclose(output_file);
}

void deleteFunction(FUNCTION *func)
// libera memoria utilizada pela struct FUNCTION
{
    free(func->initial_aps);
    _deleteResult(func->n_i);
    _deleteResult(func->n_p);
    free(func);
}
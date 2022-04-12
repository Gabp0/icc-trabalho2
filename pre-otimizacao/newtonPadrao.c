// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacoes das funcoes da biblioteca newtonPadrao.h

#include "newtonPadrao.h"
#include "utils.h"
#include "functions.h"
#include "gaussianElimination.h"
#include "Rosenbrock.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matheval.h>

NEWTON_P *initNewtonP(FUNCTION *func)
// aloca memoria para a struct NEWTON_P
{
    NEWTON_P *new = malloc(sizeof(NEWTON_P));
    if (!new)
        exitStatus(MEM_ALOC);

    // new->gradiente = malloc(sizeof(void *) * func->var_num);

    // new->hessiana = malloc(sizeof(void **) * func->var_num);
    // for (int i = 0; i < func->var_num; i++)
    //     new->hessiana[i] = malloc(sizeof(void *) * func->var_num);

    new->syst = initLS(func->var_num);
    new->X_i = copyDoubleArray(func->initial_aps, func->var_num);
    new->n = func->var_num;
    new->aprox_newtonP = calloc(sizeof(double), func->it_num + 1);
    if (!new->aprox_newtonP)
        exitStatus(MEM_ALOC);

    return new;
}

void _deleteNewtonP(NEWTON_P *np)
// libera memoria utilizada pela struct _np_
{
    free(np->aprox_newtonP);

    // for (int i = 0; i < np->n; i++)
    //     evaluator_destroy(np->gradiente[i]);

    // free(np->gradiente);

    deleteLS(np->syst);

    // for (int i = 0; i < np->n; i++)
    // {
    //     for (int j = 0; j < np->n; j++)
    //         evaluator_destroy(np->hessiana[i][j]);
    //     free(np->hessiana[i]);
    // }
    // free(np->hessiana);

    free(np->X_i);
    free(np);
}

void NewtonPadrao(FUNCTION *func)
// encontra as raizes da funcao _func_ utilizando o metodo de newton padrao
{
    func->n_p->timeFull -= timestamp();

    NEWTON_P *np = initNewtonP(func);

    func->n_p->timeDer -= timestamp();
    // Gradiente(func, np->gradiente);              // gera as funcoes do vetor gradiente
    // Hessiana(func, np->gradiente, np->hessiana); // gera as funcoes da matriz hessiana
    func->n_p->timeDer += timestamp();

    for (int k = 0; k <= func->it_num; k++) // testa numero de iteracoes
    {
        // np->aprox_newtonP[k] = evaluator_evaluate(func->evaluator, func->var_num, func->names, np->X_i); // f(X_i)
        np->aprox_newtonP[k] = rosenbrock(np->X_i, func->var_num); // f(X_i)

        for (int i = 0; i < func->var_num; i++) // gradiente f(X_i)
            // np->syst->b[i] = evaluator_evaluate(np->gradiente[i], func->var_num, func->names, np->X_i) * -1; // oposto resultado do gradiente para o calculo do sistema linear
            np->syst->b[i] = rosenbrock_dx(i, np->X_i, func->var_num) * -1; // oposto resultado do gradiente para o calculo do sistema linear

        func->n_p->it_num++; // numero de iteracoes utilizadas no metodo

        if (norma(np->syst->b, func->var_num) < func->t_ep) // testa || gradiente de f(X_i) || < eps
            break;

        for (int i = 0; i < func->var_num; i++)
            for (int j = 0; j < func->var_num; j++) // calcula a hessiana de X_i
                // np->syst->A[i][j] = evaluator_evaluate(np->hessiana[i][j], func->var_num, func->names, np->X_i);
                np->syst->A[i][j] = rosenbrock_dxdy(i, j, np->X_i, func->var_num);

        func->n_p->timeSL -= timestamp();
        gaussianElimination(np->syst); // resolve o sistema linear utilizando a eliminacao de guass
        func->n_p->timeSL += timestamp();

        for (int i = 0; i < func->var_num; i++)
            np->X_i[i] += np->syst->X[i]; // calcula X_i+1

        if (norma(np->X_i, func->var_num) < __DBL_EPSILON__) // testa || delta_i || < eps2
        {
            func->n_p->it_num++;
            break;
        }
    }

    func->n_p->f_k = copyDoubleArray(np->aprox_newtonP, func->n_p->it_num); // resultados do sistema linear
    _deleteNewtonP(np);
    func->n_p->timeFull += timestamp();
}

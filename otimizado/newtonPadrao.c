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
#include <likwid.h>

NEWTON_P *initNewtonP(FUNCTION *restrict func)
// aloca memoria para a struct NEWTON_P
{
    NEWTON_P *new = malloc(sizeof(NEWTON_P));
    if (!new)
        exitStatus(MEM_ALOC);

    new->syst = initLS(func->var_num);
    new->X_i = copyDoubleArray(func->initial_aps, func->var_num);
    new->n = func->var_num;
    new->aprox_newtonP = malloc(sizeof(double) * pad(func->it_num + 1));
    if (!new->aprox_newtonP)
        exitStatus(MEM_ALOC);
    memset(new->aprox_newtonP, 0, func->it_num + 1);

    return new;
}

void _deleteNewtonP(NEWTON_P *restrict np)
// libera memoria utilizada pela struct _np_
{

    free(np->aprox_newtonP);

    deleteLS(np->syst);

    free(np->X_i);
    free(np);
}

void NewtonPadrao(FUNCTION *restrict func)
// encontra as raizes da funcao _func_ utilizando o metodo de newton padrao
{
    // func->n_p->timeFull -= timestamp();

    NEWTON_P *np = initNewtonP(func);

    // func->n_p->timeDer -= timestamp();
    // func->n_p->timeDer += timestamp();

    string_t markerHessiana = markerName("Newton Padrao - Hessiana", func->var_num);
    string_t markerGradiente = markerName("Newton Padrao - Gradiente", func->var_num);
    string_t markerSL = markerName("NewtonPadrao - Sistema Linear", func->var_num);

    for (int k = 0; k <= func->it_num; k++) // testa numero de iteracoes
    {
        np->aprox_newtonP[k] = rosenbrock(np->X_i, func->var_num); // f(X_i)

        LIKWID_MARKER_START(markerGradiente);

        for (int i = 0; i < func->var_num; i++)                             // gradiente f(X_i)
            np->syst->b[i] = rosenbrock_dx(i, np->X_i, func->var_num) * -1; // oposto resultado do gradiente para o calculo do sistema linear

        LIKWID_MARKER_STOP(markerGradiente);

        func->n_p->it_num++; // numero de iteracoes utilizadas no metodo

        if (sq_norma(np->syst->b, func->var_num) < (func->t_ep * func->t_ep)) // testa || gradiente de f(X_i) || < eps
            break;

        LIKWID_MARKER_START(markerHessiana);
        for (int i = 0; i < func->var_num; i++)
            for (int j = 0; j < func->var_num; j++) // calcula a hessiana de X_i
                np->syst->A[i][j] = rosenbrock_dxdy(i, j, np->X_i, func->var_num);
        LIKWID_MARKER_STOP(markerHessiana);

        LIKWID_MARKER_START(markerSL);
        gaussianElimination(np->syst); // resolve o sistema linear utilizando a eliminacao de gauss
        LIKWID_MARKER_STOP(markerSL);

        double x[4];
        for (int i = 0; i < func->var_num % 4; i++)
            np->X_i[i] += np->syst->X[i];
        for (int i = func->var_num % 4; i < func->var_num; i += 4) // vetorizacao do loop
        {
            x[0] = np->X_i[i] + np->syst->X[i]; // calcula X_i+1
            x[1] = np->X_i[i + 1] + np->syst->X[i + 1];
            x[2] = np->X_i[i + 2] + np->syst->X[i + 2];
            x[3] = np->X_i[i + 3] + np->syst->X[i + 3];
            np->X_i[i] = x[0];
            np->X_i[i + 1] = x[1];
            np->X_i[i + 2] = x[2];
            np->X_i[i + 3] = x[3];
        }

        if (sq_norma(np->X_i, func->var_num) < (__DBL_EPSILON__ * __DBL_EPSILON__)) // testa || delta_i || < eps2
        {
            func->n_p->it_num++;
            break;
        }
    }

    func->n_p->f_k = copyDoubleArray(np->aprox_newtonP, func->n_p->it_num); // resultados do sistema linear
    _deleteNewtonP(np);
    // func->n_p->timeFull += timestamp();
}

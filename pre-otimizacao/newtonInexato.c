// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacao das funcoes da biblioteca newtonInexato.h

#include "newtonInexato.h"
#include "utils.h"
#include "functions.h"
#include "gaussSeidel.h"
#include "Rosenbrock.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <likwid.h>

NEWTON_I *_initNewtonI(FUNCTION *func)
// aloca memoria para a struct usada pelo metodo de newton inexato
{
    NEWTON_I *new = malloc(sizeof(NEWTON_I));
    if (!new)
        exitStatus(MEM_ALOC);

    new->syst = initLSGS(func->var_num);
    new->X_i = copyDoubleArray(func->initial_aps, func->var_num);
    new->n = func->var_num;
    new->aprox_newtonI = calloc(sizeof(double), func->it_num + 1);
    if (!new->aprox_newtonI)
        exitStatus(MEM_ALOC);

    return new;
}

void _deleteNewtonI(NEWTON_I *ni)
// libera memoria utilizada pela struct do newton inexato
{
    free(ni->aprox_newtonI);
    deleteLSGS(ni->syst);

    free(ni->X_i);
    free(ni);
}

void NewtonInexato(FUNCTION *func)
// encontra as raizes da funcao utilizando o metodo de newton inexato
{
    NEWTON_I *ni = _initNewtonI(func);
    string_t markerHessiana = markerName("Newton Inexato - Hessiana", func->var_num);
    string_t markerGradiente = markerName("Newton Inexato - Gradiente", func->var_num);
    string_t markerSL = markerName("Newton Inexato - Sistema Linear", func->var_num);

    for (int k = 0; k <= func->it_num; k++) // testa numero de iteracoes
    {
        ni->aprox_newtonI[k] = rosenbrock(ni->X_i, func->var_num); // f(X_i)

        LIKWID_MARKER_START(markerGradiente);
        for (int i = 0; i < func->var_num; i++)                             // gradiente f(X_i)
            ni->syst->b[i] = rosenbrock_dx(i, ni->X_i, func->var_num) * -1; // oposto resultado do gradiente para o calculo do sistema linear
        LIKWID_MARKER_STOP(markerGradiente);

        func->n_i->it_num++; // numero de iteracoes utilizadas no metodo

        if (norma(ni->syst->b, func->var_num) < func->t_ep) // testa || gradiente de f(X_i) || < eps
            break;

        LIKWID_MARKER_START(markerHessiana);
        for (int i = 0; i < func->var_num; i++) // calcula a hessiana de X_i
            for (int j = 0; j < func->var_num; j++)
                ni->syst->A[i][j] = rosenbrock_dxdy(i, j, ni->X_i, func->var_num);
        LIKWID_MARKER_STOP(markerHessiana);

        for (int i = 0; i < func->var_num; i++)
            ni->syst->X[i] = 0;

        LIKWID_MARKER_START(markerSL);
        gaussSeidel(ni->syst); // resolve o sistema linear utilizando gauss-seidel
        LIKWID_MARKER_STOP(markerSL);

        for (int i = 0; i < func->var_num; i++)
            ni->X_i[i] += ni->syst->X[i]; // calcula X_i+1

        if (norma(ni->X_i, func->var_num) < __DBL_EPSILON__) // testa || delta_i || < eps2
        {
            func->n_i->it_num++;
            break;
        }
    }

    func->n_i->f_k = copyDoubleArray(ni->aprox_newtonI, func->n_i->it_num); // resultado do sistema
    _deleteNewtonI(ni);
}

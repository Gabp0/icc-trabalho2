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

    // new->gradiente = malloc(sizeof(void *) * func->var_num);

    // new->hessiana = malloc(sizeof(void **) * func->var_num);
    // for (int i = 0; i < func->var_num; i++)
    //     new->hessiana[i] = malloc(sizeof(void *) * func->var_num);

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

    // for (int i = 0; i < ni->n; i++)
    // evaluator_destroy(ni->gradiente[i]);

    // free(ni->gradiente);

    deleteLSGS(ni->syst);

    // for (int i = 0; i < ni->n; i++)
    // {
    //     for (int j = 0; j < ni->n; j++)
    //         evaluator_destroy(ni->hessiana[i][j]);
    //     free(ni->hessiana[i]);
    // }
    // free(ni->hessiana);

    free(ni->X_i);
    free(ni);
}

void NewtonInexato(FUNCTION *func)
// encontra as raizes da funcao utilizando o metodo de newton inexato
{
    func->n_i->timeFull -= timestamp();

    NEWTON_I *ni = _initNewtonI(func);

    func->n_i->timeDer -= timestamp();
    // Gradiente(func, ni->gradiente);              // gera as funcoes do vetor gradiente
    // Hessiana(func, ni->gradiente, ni->hessiana); // gera as funcoes da matriz hessiana
    func->n_i->timeDer += timestamp();

    string_t *markerHessiana = markerName("NewtonInexato_Hessiana",func->var_num);
    string_t *markerGradiente = markerName("NewtonInexato_Gradiente",func->var_num);
    string_t *markerSL = markerName("NewtonInexato_SL",func->var_num);

    for (int k = 0; k <= func->it_num; k++) // testa numero de iteracoes
    {
        // ni->aprox_newtonI[k] = evaluator_evaluate(func->evaluator, func->var_num, func->names, ni->X_i); // f(X_i)
        ni->aprox_newtonI[k] = rosenbrock(ni->X_i, func->var_num); // f(X_i)

        LIKWID_MARKER_START(markerGradiente);
        for (int i = 0; i < func->var_num; i++) // gradiente f(X_i)
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

        func->n_i->timeSL -= timestamp();
        LIKWID_MARKER_START(markerSL);
        gaussSeidel(ni->syst); // resolve o sistema linear utilizando gauss-seidel
        LIKWID_MARKER_STOP(markerSL);
        func->n_i->timeSL += timestamp();

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
    func->n_i->timeFull += timestamp();
}

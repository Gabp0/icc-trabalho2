// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacoes das funcoes da biblioteca newtonModificado.h

#include "utils.h"
#include "functions.h"
#include "luDecomposition.h"
#include "newtonModificado.h"
#include "Rosenbrock.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

NEWTON_M *_initNewtonM(FUNCTION *func)
// aloca memoria para a struct NEWTOM_M
{
    NEWTON_M *new = malloc(sizeof(NEWTON_M));

    if (!new)
        exitStatus(MEM_ALOC);

    // new->gradiente = malloc(sizeof(void *) * func->var_num);

    // new->hessiana = malloc(sizeof(void **) * func->var_num);
    // for (int i = 0; i < func->var_num; i++)
    //     new->hessiana[i] = malloc(sizeof(void *) * func->var_num);

    new->syst = initLSLU(func->var_num);
    new->X_i = copyDoubleArray(func->initial_aps, func->var_num);
    new->n = func->var_num;
    new->aprox_newtonM = calloc(sizeof(double), func->it_num + 1);
    if (!new->aprox_newtonM)
        exitStatus(MEM_ALOC);

    return new;
}

void _deleteNewtonM(NEWTON_M *nm)
// libera memoria utilizada pela struct _nm_
{
    free(nm->aprox_newtonM);

    // for (int i = 0; i < nm->n; i++)
    //     evaluator_destroy(nm->gradiente[i]);
    // free(nm->gradiente);

    deleteLSLU(nm->syst);

    // for (int i = 0; i < nm->n; i++)
    // {
    //     for (int j = 0; j < nm->n; j++)
    //         evaluator_destroy(nm->hessiana[i][j]);
    //     free(nm->hessiana[i]);
    // }
    // free(nm->hessiana);

    free(nm->X_i);
    free(nm);
}

void NewtonModificado(FUNCTION *func)
// encontra as raizes da funcao _func_ utilizando o metodo de newon modificado
{
    func->n_m->timeFull -= timestamp();
    double soma = 0;
    NEWTON_M *nm = _initNewtonM(func);

    func->n_m->timeDer -= timestamp();
    // Gradiente(func, nm->gradiente);              // gera as funcoes do vetor gradiente
    // Hessiana(func, nm->gradiente, nm->hessiana); // gera as funcoes da matriz hessiana
    func->n_m->timeDer += timestamp();

    soma += pow(nm->syst->X[0], 2);

    for (int k = 0; k <= func->it_num; k++) // testa numero de iteracoes
    {
        // nm->aprox_newtonM[k] = evaluator_evaluate(func->evaluator, func->var_num, func->names, nm->X_i); // f(X_i)
        nm->aprox_newtonM[k] = rosenbrock(nm->X_i, func->var_num); // f(X_i)

        for (int i = 0; i < func->var_num; i++) // gradiente f(X_i)
            // nm->syst->b[i] = evaluator_evaluate(nm->gradiente[i], func->var_num, func->names, nm->X_i) * -1; // oposto resultado do gradiente para o calculo do sistema linear
            nm->syst->b[i] = rosenbrock_dx(i, nm->X_i, func->var_num) * -1; // oposto resultado do gradiente para o calculo do sistema linear

        func->n_m->it_num++; // numero de iteracoes utilizadas no metodo

        if (norma(nm->syst->b, func->var_num) < func->t_ep) // testa || gradiente de f(X_i) || < eps
            break;

        if ((k % func->var_num) == 0)
        {
            for (int i = 0; i < func->var_num; i++)
                for (int j = 0; j < func->var_num; j++) // calcula a hessiana de X_i apenas a cada HESS_STEPS
                    // nm->syst->U[i][j] = evaluator_evaluate(nm->hessiana[i][j], func->var_num, func->names, nm->X_i);
                    nm->syst->U[i][j] = rosenbrock_dxdy(i, j, nm->X_i, func->var_num);

            factorize(nm->syst); // separa a matriz H_x no sistema linear com LU
        }

        func->n_m->timeSL -= timestamp();
        solveLU(nm->syst); // resolve o sistema linear utilizando fatoracao LU
        func->n_m->timeSL += timestamp();

        soma = 0;
        for (int i = 0; i < func->var_num; i++)
            nm->X_i[i] += nm->syst->X[i]; // calcula X_i+1

        if (norma(nm->X_i, func->var_num) < __DBL_EPSILON__) // testa || delta_i || < eps
        {
            func->n_m->it_num++;
            break;
        }
    }

    func->n_m->f_k = copyDoubleArray(nm->aprox_newtonM, func->n_m->it_num); // resultado do sistema
    _deleteNewtonM(nm);
    func->n_m->timeFull += timestamp();
}

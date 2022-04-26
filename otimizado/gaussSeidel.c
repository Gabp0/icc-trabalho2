// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacao das funcoes da biblioteca gaussSeidel.h

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matheval.h>
#include "gaussSeidel.h"

LINEAR_SYST_GS *initLSGS(int size)
// Aloca memoria para o sistema linear de tamanho size
{
    LINEAR_SYST_GS *new = malloc(sizeof(LINEAR_SYST_GS));
    if (!new)
        exitStatus(MEM_ALOC);

    new->size = size;
    int p_size = pad(size); // padding no tamanho para evitar cache trashing

    new->A = malloc(sizeof(double) * p_size * p_size);
    if (!new->A)
        exitStatus(MEM_ALOC);

    new->b = malloc(sizeof(double) * p_size);
    if (!new->b)
        exitStatus(MEM_ALOC);

    new->X = malloc(sizeof(double) * p_size);
    if (!new->X)
        exitStatus(MEM_ALOC);

    new->Xk_m1 = malloc(sizeof(double) * p_size);
    if (!new->Xk_m1)
        exitStatus(MEM_ALOC);

    return new;
}

void gaussSeidel(LINEAR_SYST_GS *restrict syst)
// resolve o sistema linear utilizando o metodo de Gauss-Seidel
{
    double soma;
    int i_m;

    for (int i = 0; i < syst->size; i++) // primeira iteracao fora do loop para remocao do if
    {
        i_m = i * syst->size;
        soma = 0;
        for (int j = 0; j < i; j++) // quebra do loop em dois eliminando o if
            soma += syst->A[i_m + j] * syst->X[j];
        for (int j = i + 1; j < syst->size; j++)
            soma += syst->A[i_m + j] * syst->X[j];

        syst->Xk_m1[i] = syst->X[i]; // guarda x[k - 1]
        syst->X[i] = (syst->b[i] - soma) / syst->A[i_m + i];
    }

    for (int k = 1; (k < IT_MAX) && (fabs(norma(syst->X, syst->size) - norma(syst->Xk_m1, syst->size)) < TOL); k++) // numero de iteracoes
    {

        for (int i = 0; i < syst->size; i++)
        {
            i_m = i * syst->size;
            soma = 0;
            for (int j = 0; j < i; j++) // quebra do loop em dois eliminando o if
                soma += syst->A[i_m + j] * syst->X[j];
            for (int j = i + 1; j < syst->size; j++)
                soma += syst->A[i_m + j] * syst->X[j];

            syst->Xk_m1[i] = syst->X[i]; // guarda x[k - 1]
            syst->X[i] = (syst->b[i] - soma) / syst->A[i_m + i];
        }
    }

    return;
}

void deleteLSGS(LINEAR_SYST_GS *restrict syst)
// libera memoria utilizada pelo sistema linear _syst_
{
    if (!syst)
        exitStatus(INV_POINTER);

    free(syst->A);
    free(syst->b);
    free(syst->X);
    free(syst->Xk_m1);
    free(syst);
}

// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacao das funcoes da biblioteca gaussSeidel.h

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matheval.h>
#include <string.h>
#include "gaussSeidel.h"

#define NDSIMD 4

LINEAR_SYST_GS *initLSGS(int size)
// Aloca memoria para o sistema linear de tamanho size
{
    LINEAR_SYST_GS *new = malloc(sizeof(LINEAR_SYST_GS));
    if (!new)
        exitStatus(MEM_ALOC);

    new->size = size;
    int p_size = pad(size); // padding no tamanho para evitar cache trashing

    new->A = malloc(sizeof(double) * size * p_size);
    if (!new->A)
        exitStatus(MEM_ALOC);

    new->b = malloc(sizeof(double) * p_size);
    if (!new->b)
        exitStatus(MEM_ALOC);

    new->X = malloc(sizeof(double) * p_size);
    if (!new->X)
        exitStatus(MEM_ALOC);
    memset(new->X, 0, sizeof(double) * p_size);

    new->Xk_m1 = malloc(sizeof(double) * p_size);
    if (!new->Xk_m1)
        exitStatus(MEM_ALOC);

    return new;
}

void gaussSeidel(LINEAR_SYST_GS *restrict syst)
// resolve o sistema linear utilizando o metodo de Gauss-Seidel
{
    double so[NDSIMD];
    double soma;
    int i_m;

    for (int k = 0; k < IT_MAX && (k == 0 || fabs(sq_norma(syst->X, syst->size) - sq_norma(syst->Xk_m1, syst->size)) > (TOL * TOL)); k++) // numero de iteracoes
    {
        for (int i = 0; i < syst->size; i++)
        {
            i_m = i * syst->size;
            soma = 0.0;
            so[0] = 0;
            so[1] = 0;
            so[2] = 0;
            so[3] = 0;
            int j;

            for (j = 0; j < i % NDSIMD; j++) // residuo
                soma += syst->A[i_m + j] * syst->X[j];

            for (; j < i; j += NDSIMD) // vetorizacao
            {
                so[0] += syst->A[i_m + j] * syst->X[j];
                so[1] += syst->A[i_m + j + 1] * syst->X[j + 1];
                so[2] += syst->A[i_m + j + 2] * syst->X[j + 2];
                so[3] += syst->A[i_m + j + 3] * syst->X[j + 3];
            }
            // quebra do loop em dois eliminando o if
            for (j = i + 1; j < i + 1 + (syst->size - i - 1) % NDSIMD; j++) // residuo
                soma += syst->A[i_m + j] * syst->X[j];

            for (; j < syst->size; j += NDSIMD) // vetorizacao
            {
                so[0] += syst->A[i_m + j] * syst->X[j];
                so[1] += syst->A[i_m + j + 1] * syst->X[j + 1];
                so[2] += syst->A[i_m + j + 2] * syst->X[j + 2];
                so[3] += syst->A[i_m + j + 3] * syst->X[j + 3];
            }

            soma += so[0] + so[1] + so[2] + so[3];

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

// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacao das funcoes da biblioteca luDecomposition.h

#include "luDecomposition.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "utils.h"

LINEAR_SYST_LU *initLSLU(int size)
// Aloca memoria para o sistema linear LU de tamanho size
{
    LINEAR_SYST_LU *new = malloc(sizeof(LINEAR_SYST_LU));
    if (!new)
        exitStatus(MEM_ALOC);

    new->size = size;

    new->L = malloc(sizeof(double *) * size);
    if (!new->L)
        exitStatus(MEM_ALOC);
    for (int i = 0; i < size; i++)
    {
        new->L[i] = malloc(sizeof(double) * size);
        if (!new->L[i])
            exitStatus(MEM_ALOC);
    }

    new->U = malloc(sizeof(double *) * size);
    if (!new->U)
        exitStatus(MEM_ALOC);
    for (int i = 0; i < size; i++)
    {
        new->U[i] = malloc(sizeof(double) * size);
        if (!new->U[i])
            exitStatus(MEM_ALOC);
    }

    new->b = malloc(sizeof(double) * size);
    if (!new->b)
        exitStatus(MEM_ALOC);

    new->X = calloc(sizeof(double), size);
    if (!new->X)
        exitStatus(MEM_ALOC);

    new->Y = malloc(sizeof(double) * size);
    if (!new->Y)
        exitStatus(MEM_ALOC);

    new->swaps = 0;
    new->swap_index = malloc(sizeof(INDEXES) * size);

    return new;
}

void _swapBLines(LINEAR_SYST_LU *syst)
// realiza as trocas feitas na matriz de coeficientes no vetor resultado
{
    double aux;

    for (int i = 0; i < syst->swaps; i++)
    {
        aux = syst->b[syst->swap_index[i].ia];
        syst->b[syst->swap_index[i].ia] = syst->b[syst->swap_index[i].ib];
        syst->b[syst->swap_index[i].ib] = aux;
    }
}

void _pivotLU(LINEAR_SYST_LU *syst, int i)
// faz o pivoteamento parcial do sistema linear e guarda as trocas feitas para o vetor resultado
{
    double max = fabs(syst->U[i][i]);
    int max_i = i;
    for (int j = i + 1; j < syst->size; ++j) // encontra o maior valor
    {
        double v = fabs(syst->U[j][i]);
        if (v > max)
        {
            max = v;
            max_i = j;
        }
    }

    if (max_i != i) // faz a substituicao
    {
        double *tmp = syst->U[i];
        syst->U[i] = syst->U[max_i];
        syst->U[max_i] = tmp;

        syst->swap_index[syst->swaps].ia = i; // guarda os indices da substituicao
        syst->swap_index[syst->swaps].ib = max_i;
        syst->swaps++;
    }
}

void factorize(LINEAR_SYST_LU *syst)
// divide a matriz de coeficientes na triangular superior U e triangular inferior L
{
    syst->swaps = 0;
    for (int i = 0; i < syst->size; ++i)
    {
        _pivotLU(syst, i); // primeira linha com o maior valor

        syst->L[i][i] = 1.0; // diagonal principal da L

        for (int k = i + 1; k < syst->size; ++k)
        {
            syst->L[k][i] = syst->U[k][i] / syst->U[i][i]; // U recebe os coeficientes da triangulacao

            if (!isValidNum(syst->L[k][i]))
                exitStatus(ZERO_DIV);

            syst->U[k][i] = 0.0;

            for (int j = i + 1; j < syst->size; ++j)
                syst->U[k][j] -= syst->U[i][j] * syst->L[k][i];
        }
    }
}

void _subsLU(LINEAR_SYST_LU *syst)
// encontra o valor de Y na fatoracao LU substituindo na matriz L a partir da primeira linha
{
    _swapBLines(syst);
    for (int i = 0; i < syst->size; ++i)
    {
        syst->Y[i] = syst->b[i];
        for (int j = 0; j < i; j++)
            syst->Y[i] -= syst->L[i][j] * syst->Y[j];
    }
}

void _retrossubsLU(LINEAR_SYST_LU *syst)
// encontra o valor de X na fatoracao LU substituindo na matriz U a partir da utima linha
{
    for (int i = syst->size - 1; i >= 0; --i)
    {
        syst->X[i] = syst->Y[i];
        for (int j = i + 1; j < syst->size; j++)
            syst->X[i] -= syst->U[i][j] * syst->X[j];
        syst->X[i] /= syst->U[i][i];
    }
}

void solveLU(LINEAR_SYST_LU *syst)
{
    _subsLU(syst);
    _retrossubsLU(syst);
}

void deleteLSLU(LINEAR_SYST_LU *syst)
// libera a memoria utilizada pelo sistema linear LU _syst_
{
    if (!syst)
        exitStatus(INV_POINTER);

    for (int i = 0; i < syst->size; i++)
        free(syst->L[i]);
    free(syst->L);
    for (int i = 0; i < syst->size; i++)
        free(syst->U[i]);
    free(syst->U);
    free(syst->b);
    free(syst->X);
    free(syst->Y);
    free(syst->swap_index);
    free(syst);
}

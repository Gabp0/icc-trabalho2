// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacoes das funcoes da biblioteca gaussianElimination.h

#include "gaussianElimination.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

LINEAR_SYST *initLS(int size)
// Aloca memoria para o sistema linear de tamanho size
{
  LINEAR_SYST *new = malloc(sizeof(LINEAR_SYST));
  if (!new)
    exitStatus(MEM_ALOC);

  new->size = size;
  new->A = initDoubleMatrix(size);

  new->b = malloc(sizeof(double) * pad(size));
  if (!new->b)
    exitStatus(MEM_ALOC);

  new->X = malloc(sizeof(double) * pad(size));
  if (!new->X)
    exitStatus(MEM_ALOC);

  return new;
}

void deleteLS(LINEAR_SYST *restrict syst)
// Libera a memoria utilizada pelo sistema linear syst
{
  if (!syst)
    exitStatus(INV_POINTER);

  free(syst->A);
  free(syst->b);
  free(syst->X);
  free(syst);
}

void _pivot(LINEAR_SYST *restrict syst, int i)
// faz o pivoteamento parcial do sistema linear _syst_
{
  double max = fabs(syst->A[i][i]);
  int max_i = i;

  for (int j = i + 1; j < syst->size; ++j) // encontra o maior valor da coluna
  {
    double v = fabs(syst->A[j][i]);
    if (v > max)
    {
      max = v;
      max_i = j;
    }
  }

  if (max_i != i) // substitui
  {
    double *tmp = syst->A[i];
    syst->A[i] = syst->A[max_i];
    syst->A[max_i] = tmp;

    double aux = syst->b[i];
    syst->b[i] = syst->b[max_i];
    syst->b[max_i] = aux;
  }
}

void _retrossubs(LINEAR_SYST *restrict syst)
// encontra os valores de X substituindo a partir da ultima linha do sl
{
  // unroll & jam do loop
  // for (int i = syst->size - 1; i >= 0; --i)
  //    {
  //      syst->X[i] = syst->b[i];
  //      for (int j = i + 1; j < syst->size; j++)
  //        syst->X[i] -= syst->A[i][j] * syst->X[j];
  //      syst->X[i] /= syst->A[i][i];

  int aux = (syst->size - (syst->size % 4));

  for (int i = (syst->size - 1); i > aux - 1; --i)
  {
    syst->X[i] = syst->b[i];
    for (int j = i + 1; j < syst->size; j++)
      syst->X[i] -= syst->A[i][j] * syst->X[j];
    syst->X[i] /= syst->A[i][i];
  }

  for (int i = aux - 1; i >= 0; i -= 4)
  {
    syst->X[i] = syst->b[i];
    syst->X[i - 1] = syst->b[i - 1];
    syst->X[i - 2] = syst->b[i - 2];
    syst->X[i - 3] = syst->b[i - 3];
    for (int j = i + 1; j < syst->size; j++)
    {
      syst->X[i] -= syst->A[i][j] * syst->X[j];
      syst->X[i - 1] -= syst->A[i - 1][j] * syst->X[j];
      syst->X[i - 2] -= syst->A[i - 2][j] * syst->X[j];
      syst->X[i - 3] -= syst->A[i - 3][j] * syst->X[j];
    }
    syst->X[i] /= syst->A[i][i];
    syst->X[i - 1] /= syst->A[i - 1][i - 1];
    syst->X[i - 2] /= syst->A[i - 2][i - 2];
    syst->X[i - 3] /= syst->A[i - 3][i - 3];
  }
}

void _triang(LINEAR_SYST *restrict syst)
// coloca o sistema na forma escada
{
  for (int i = 0; i < syst->size; ++i)
  {
    _pivot(syst, i); // primeira linha tem o maior valor

    for (int k = i + 1; k < syst->size; ++k)
    {
      double m = syst->A[k][i] / syst->A[i][i];

      // if (isnan(m) || isinf(m))
      //   exitStatus(ZERO_DIV);

      syst->A[k][i] = 0.0;
      for (int j = i + 1; j < syst->size; ++j)
        syst->A[k][j] -= syst->A[i][j] * m;

      syst->b[k] -= syst->b[i] * m;
    }
  }
}

void gaussianElimination(LINEAR_SYST *restrict syst)
// resolve o sistema linear utilizando a eliminacao de gauss
{
  _triang(syst);
  _retrossubs(syst);
}

// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacoes das funcoes da biblioteca gaussianElimination.h

#include "gaussianElimination.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"

LINEAR_SYST *initLS(int size)
// Aloca memoria para o sistema linear de tamanho size
{
	LINEAR_SYST *new = malloc(sizeof(LINEAR_SYST));
	if (!new)
		exitStatus(MEM_ALOC);

	new->size = size;
	int p_size = pad(size); // padding no tamanho para evitar cache trashing

	new->A = malloc(p_size * sizeof(double *));
	if (!new->A)
		exitStatus(MEM_ALOC);
	new->A[0] = malloc(p_size * p_size * sizeof(double)); // aloca um vetor com todos os elementos da matriz
	if (!new->A[0])
		exitStatus(MEM_ALOC);
	for (int i = 1; i < size; i++) // ajusta os demais ponteiros de linhas (i > 0)
		new->A[i] = new->A[0] + i *size;

	new->b = malloc(sizeof(double) * p_size);
	if (!new->b)
		exitStatus(MEM_ALOC);

	new->X = malloc(sizeof(double) * p_size);
	if (!new->X)
		exitStatus(MEM_ALOC);

	return new;
}

void deleteLS(LINEAR_SYST *restrict syst)
// Libera a memoria utilizada pelo sistema linear syst
{
	if (!syst)
		exitStatus(INV_POINTER);

	free(syst->A[0]);
	free(syst->A);

	free(syst->b);
	// printf("ok\n");
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

#define NDSIMD 4

void _retrossubs(LINEAR_SYST *restrict syst)
// encontra os valores de X substituindo a partir da ultima linha do sl
{
	double x[4];
	int j;

	// loop original
	// for (int i = syst->size - 1; i >= 0; --i)
	// {
	// 	syst->X[i] = syst->b[i];
	// 	for (int j = i + 1; j < syst->size; j++)
	// 		syst->X[i] -= syst->A[i][j] * syst->X[j];
	// 	syst->X[i] /= syst->A[i][i];
	// }

	memcpy(syst->X, syst->b, sizeof(double) * syst->size);
	for (int i = syst->size - 1; i >= 0; i--)
	{
		x[0] = 0.0;
		x[1] = 0.0;
		x[2] = 0.0;
		x[3] = 0.0;

		for (j = i + 1; j < (syst->size % NDSIMD); j++)
			x[0] += syst->A[i][j] * syst->X[j];

		for (; j < syst->size; j += NDSIMD) // unroll & jam do loop
		{
			x[0] += syst->A[i][j] * syst->X[j];
			x[1] += syst->A[i][j + 1] * syst->X[j + 1];
			x[2] += syst->A[i][j + 2] * syst->X[j + 2];
			x[3] += syst->A[i][j + 3] * syst->X[j + 3];
		}
		syst->X[i] -= (x[0] + x[1] + x[2] + x[3]);
		syst->X[i] /= syst->A[i][i];
	}
}

void _triang(LINEAR_SYST *restrict syst)
// coloca o sistema na forma escada
{
	double m;

	for (int i = 0; i < syst->size; ++i)
	{
		_pivot(syst, i); // primeira linha tem o maior valor

		for (int k = i + 1; k < syst->size; ++k)
		{
			m = syst->A[k][i] / syst->A[i][i];

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

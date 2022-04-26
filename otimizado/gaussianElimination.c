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

void _retrossubs(LINEAR_SYST *restrict syst)
// encontra os valores de X substituindo a partir da ultima linha do sl
{
	// unroll & jam do loop
	for (int i = syst->size - 1; i >= 0; --i)
	{
		syst->X[i] = syst->b[i];
		for (int j = i + 1; j < syst->size; j++)
			syst->X[i] -= syst->A[i][j] * syst->X[j];
		syst->X[i] /= syst->A[i][i];
	}

	// int size_strd = (syst->size - (syst->size % 4));

	// for (int i = (syst->size - 1); i > size_strd - 1; --i) // residuo
	// {
	// 	syst->X[i] = syst->b[i];
	// 	for (int j = i + 1; j < syst->size; j++)
	// 		syst->X[i] -= syst->A[i][j] * syst->X[j];
	// 	syst->X[i] /= syst->A[i][i];
	// }

	// for (int i = size_strd - 1; i >= 0; i -= 4) // unroll e jam com 4 de stride
	// {
	// 	syst->X[i] = syst->b[i];
	// 	syst->X[i - 1] = syst->b[i - 1];
	// 	syst->X[i - 2] = syst->b[i - 2];
	// 	syst->X[i - 3] = syst->b[i - 3];
	// 	for (int j = i + 1; j < syst->size; j++)
	// 	{
	// 		syst->X[i] -= syst->A[i][j] * syst->X[j];
	// 		syst->X[i - 1] -= syst->A[i - 1][j] * syst->X[j];
	// 		syst->X[i - 2] -= syst->A[i - 2][j] * syst->X[j];
	// 		syst->X[i - 3] -= syst->A[i - 3][j] * syst->X[j];
	// 	}
	// 	syst->X[i] /= syst->A[i][i];
	// 	syst->X[i - 1] /= syst->A[i - 1][i - 1];
	// 	syst->X[i - 2] /= syst->A[i - 2][i - 2];
	// 	syst->X[i - 3] /= syst->A[i - 3][i - 3];
	// }
}

void _triang(LINEAR_SYST *restrict syst)
// coloca o sistema na forma escada
{
	double m;
	// int size_strd = syst->size - (syst->size % 4);

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

		// for (int k = i + 1; k < size_strd - 1; k += 4)
		// {
		// 	// printf("ok\n");

		// 	m[0] = syst->A[k][i] / syst->A[i][i];
		// 	syst->A[k][i] = 0.0;
		// 	m[1] = syst->A[k + 1][i] / syst->A[i][i];
		// 	syst->A[k + 1][i] = 0.0;
		// 	m[2] = syst->A[k + 2][i] / syst->A[i][i];
		// 	syst->A[k + 2][i] = 0.0;
		// 	// printf("aqui %d \n", k);
		// 	m[3] = syst->A[k + 3][i] / syst->A[i][i];
		// 	syst->A[k + 3][i] = 0.0;

		// 	// printf("%d %d\n", k, i);

		// 	for (int j = i + 1; j < syst->size; ++j)
		// 	{
		// 		syst->A[k][j] -= syst->A[i][j] * m[0];
		// 		syst->A[k + 1][j] -= syst->A[i][j] * m[1];
		// 		syst->A[k + 2][j] -= syst->A[i][j] * m[2];
		// 		syst->A[k + 3][j] -= syst->A[i][j] * m[3];
		// 	}

		// 	syst->b[k] -= syst->b[i] * m[0];
		// 	syst->b[k + 1] -= syst->b[i] * m[1];
		// 	syst->b[k + 2] -= syst->b[i] * m[2];
		// 	syst->b[k + 3] -= syst->b[i] * m[3];
		// }

		// for (int k = max(i + 1, size_strd - 1, 0); k < syst->size; k += 1) // residuo
		// {

		// 	m[0] = syst->A[k][i] / syst->A[i][i];
		// 	syst->A[k][i] = 0.0;

		// 	for (int j = i + 1; j < syst->size; ++j)
		// 		syst->A[k][j] -= syst->A[i][j] * m[0];

		// 	syst->b[k] -= syst->b[i] * m[0];
		// }
	}
}

void gaussianElimination(LINEAR_SYST *restrict syst)
// resolve o sistema linear utilizando a eliminacao de gauss
{
	_triang(syst);
	_retrossubs(syst);
}

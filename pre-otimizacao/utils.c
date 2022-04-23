// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Implementacoes das funcoes da biblioteca utils.h

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include "utils.h"

void exitStatus(int code)
// encerra o programa com o determinado status de saida e mensagem de erro
{
	switch (code)
	{
	case MEM_ALOC:
		perror("Falha na alocação de memória");
		exit(code);

	case INV_POINTER:
		perror("Ponteiro inválido");
		exit(code);

	case MATHEVAL_ERR:
		perror("Erro ao gerar avaliador para a libmatheval");
		exit(code);

	case ARG_NUM:
		perror("Número de argumentos na chamada do programa inválido");
		exit(code);

	case ARG_INV:
		perror("Argumento inválido");
		exit(code);

	case FOPEN_ERR:
		perror("Não foi possível abrir o arquivo");
		exit(code);

	case ZERO_DIV:
		perror("Divisão por zero");
		exit(code);

	default:
		break;
	}
}

double timestamp(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double)(tp.tv_sec * 1000.0 + tp.tv_usec / 1000.0));
}

void prnVetorFloat(float *v, int n)
{
	for (int i = 0; i < n; ++i)
		printf("%g ", v[i]);
	printf("\n");
}

void prnVetorDouble(double *v, int n)
{
	for (int i = 0; i < n; ++i)
		fprintf(stderr, "%lg ", v[i]);
	fprintf(stderr, "\n");
}

void prnVetorLongDouble(long double *v, int n)
{
	for (int i = 0; i < n; ++i)
		printf("%Lg ", v[i]);
	printf("\n");
}

double *copyDoubleArray(double *a, int size)
// aloca um array de double, copia o conteudo de _a_ e retorna o ponteiro
{
	double *new = malloc(sizeof(double) * size);
	for (int i = 0; i < size; i++)
		new[i] = a[i];

	return new;
}

double **initDoubleMatrix(int size)
// aloca memoria para uma matriz de double de tamanho _size_ * _size_ e retorna o ponteiro para ela
{
	double **A = malloc(sizeof(double *) * size);
	if (!A)
		exitStatus(MEM_ALOC);

	for (int i = 0; i < size; i++)
	{
		A[i] = malloc(sizeof(double) * size);
		if (!A[i])
			exitStatus(MEM_ALOC);
	}

	return A;
}

char *getArgs(int argc, char **argv)
// separa os argumentos passados pela linha de comando
{
	if (argc == 3)
	{
		if (!strcmp(argv[1], "-o"))
			return argv[2];

		else
			exitStatus(ARG_INV);
	}
	else if (argc == 1)
		return NULL;

	exitStatus(ARG_NUM);
	return NULL;
}

int max(int a, int b, int c)
// maximo entre tres int
{
	int aux = a > b ? a : b;
	return aux > c ? aux : c;
}

double norma(double *array, int size)
// norma euclidiana do vetor
{
	double soma = 0;
	for (int i = 0; i < size; i++)
		soma = array[i] * array[i];
	return sqrt(soma);
}

int isValidNum(double num)
// testa se o numero nao é nan ou inf
{
	return !(isnan(num) || isinf(num));
}

string_t markerName(string_t baseName, int n)
{
  string_t mark = (string_t) malloc( (strlen(baseName)+1) + (log10(n)+1) + 1 );

  sprintf(mark, "%s_%u", baseName,n);

  // printf("*** %s\n", mark);

  return mark;

}

int isPot2(int n)
{
  int k;
  return (k = log2(n)) == log2(n) ;
}

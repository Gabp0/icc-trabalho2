// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Biblioteca gaussSeidel.h
// Funcoes para representacao e resolucao de sistemas lineares com utilizando o metodo de Gauss-Seidel

#ifndef __GAUSSSEILDEL__
#define __GAUSSSEILDEL__

#define IT_MAX 50    // numero maximo de iteracoes do metodo de Gauss-Seidel
#define TOL 0.000001 // tolerancia minima para o metodo de Gauss-Seidel

typedef struct linear_syst_gs
// struct para representacao de um sistema linear
{
    double **A;    // matriz de coeficientes A
    double *X;     // vetor das variaveis
    double *Xk_m1; // vetor das variaveis da iteracao anterior
    double *b;     // vetor do resultado
    int size;      // tamanho
} LINEAR_SYST_GS;

LINEAR_SYST_GS *initLSGS(int size);              // Aloca memoria para um sistema linear de tamanho _size_ para resolucao por Gauss-Seidel
void gaussSeidel(LINEAR_SYST_GS *restrict syst); // Resolve o sistema linear _syst_ utilizando o metodo de Gauss-Seidel
void deleteLSGS(LINEAR_SYST_GS *restrict syst);  // Libera memoria utilizada pelo sistema linear _syst_

#endif

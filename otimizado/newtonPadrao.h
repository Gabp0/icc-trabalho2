// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
//  Biblioteca newtonPadrao.h
//  Encontra as raizes da funcao utilizando o metodo de Newton Padrao com eliminacao de Gauss

#ifndef __NEWTONPADRAO__
#define __NEWTONPADRAO__

#include "functions.h"
#include "gaussianElimination.h"

typedef struct newton_padrao
{
    // void ***hessiana;      // matriz de equacoes da hessiana
    // void **gradiente;      // vetor de equacoes do gradiente
    LINEAR_SYST *syst;     // sistema linear H(x_i) * delta = - Gf(X_i)
    double *X_i;           // vetor de solucoes do polinomio
    int n;                 // numero de variaveis
    double *aprox_newtonP; // vetor f(X_k)
    double eps;            // tolerancia epsilon
} NEWTON_P;

void NewtonPadrao(FUNCTION *restrict func); // Encontra as raizes da funcao _func_ utilizando o metodo de newton padrao

#endif
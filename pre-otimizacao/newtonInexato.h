// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Biblioteca newtonInexato.h
// Encontra as raizes da funcao utilizando o metodo de Newton Inexato com Gauss-Seidel

#ifndef __NEWTONINEXATO__
#define __NEWTONINEXATO__

#include "functions.h"
#include "gaussSeidel.h"

typedef struct newton_inexato
{
    // void ***hessiana;      // matriz de equacoes da hessiana
    // void **gradiente;      // vetor de equacoes do gradiente
    LINEAR_SYST_GS *syst;  // sistema linear H(x_i) * delta = - Gf(X_i)
    double *X_i;           // vetor de solucoes do polinomio
    int n;                 // numero de variaveis
    double *aprox_newtonI; // vetor com f(x_k)
    double eps;            // tolerancia epsilon
} NEWTON_I;

void NewtonInexato(FUNCTION *func); // Encontra as raizes da funcao _func_ utilizando o metodo de newton inexato

#endif
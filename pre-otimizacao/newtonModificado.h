// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Biblioteca newtonModificado.h
// Encontra as raizes da funcao utilizando o metodo de Newton Modificado com Fatoracao LU

#ifndef __NEWTONMODIFICADO__
#define __NEWTONMODIFICADO__

#include "functions.h"
#include "luDecomposition.h"

#define HESS_STEPS 5

typedef struct newton_modificado
{
    // void ***hessiana;      // matriz de equacoes da hessiana
    // void **gradiente;      // vetor de equacoes do gradiente
    LINEAR_SYST_LU *syst;  // sistema linear H'(X_i) * delta = - Gf(X_i)
    double *X_i;           // vetor de solucoes do polinomio
    int n;                 // numero de variaveis
    double *aprox_newtonM; // vetor f(X_k)
    double eps;            // tolerancia epsilon
} NEWTON_M;

void NewtonModificado(FUNCTION *func); // Encontra as raizes da funcao _func_ utilizando o metodo de newton modificado

#endif
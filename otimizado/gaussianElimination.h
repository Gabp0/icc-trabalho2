// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
//  Biblioteca gaussianElimination.h
//  Funcoes para representacao e resolucao de sistemas lineares utilizando o metodo da eliminacao de Gauss

#ifndef __GAUSS__
#define __GAUSS__

typedef struct linear_syst
// struct para representacao de um sistema linear
{
    double **A; // matriz de coeficientes A
    double *X;  // vetor das variaveis
    double *b;  // vetor do resultado
    int size;   // tamanho
} LINEAR_SYST;

LINEAR_SYST *initLS(int size);                        // Aloca memoria para um sistema linear de tamanho _size_ e retorna um ponteiro para ele
void deleteLS(LINEAR_SYST *restrict syst);            // Libera memoria utilizada pelo sistema linear _syst_
void gaussianElimination(LINEAR_SYST *restrict syst); // Resolve o sistema linear _syst_ utilizando o metodo da eliminacao de Gauss

#endif

// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092

#ifndef __LU_DEC__
#define __LU_DEC__

typedef struct indexes
{
    int ia;
    int ib;
} INDEXES;

// struct para representacao de um sistema linear
typedef struct linear_syst_lu
{
    double **L;          // matriz de coeficientes fatorados L
    double **U;          // matriz de coeficientes fatorados U
    double *X;           // vetor das variaveis
    double *Y;           // Ly = b
    double *b;           // vetor do resultado
    INDEXES *swap_index; // vetor de index para as trocas do pivoteamento
    int swaps;           // quantidade de trocas
    int size;            // tamanho
} LINEAR_SYST_LU;

LINEAR_SYST_LU *initLSLU(int size);    // Aloca memoria para um sistema linear do tamanho _size_ para resolucao por fatoracao LU
void deleteLSLU(LINEAR_SYST_LU *syst); // Libera memoria utilizada pelo sistema linear _syst_
void factorize(LINEAR_SYST_LU *syst);  // Fatora a matriz armazenada em syst->U na triangular inferior L e triangular superior U
void solveLU(LINEAR_SYST_LU *syst);    // Resolve o sistema linear para fatoracao LU

#endif
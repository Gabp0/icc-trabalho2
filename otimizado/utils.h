// Gabriel de Oliveira Pontarolo GRR20203895
// Rodrigo Saviam Soffner GRR20205092
// Biblioteca utils.h
// Funcoes de uso geral no programa

#ifndef __UTILS__
#define __UTILS__

#include <stdlib.h>
#include <sys/time.h>

// codigos de saida de erro
#define MEM_ALOC 101     // falha de alocacao de memoria
#define INV_POINTER 102  // ponteiro null
#define MATHEVAL_ERR 103 // ponteiro null como retorno de funcao da lib matheval
#define ARG_NUM 104      // numero de argumentos passados pela linha de comando invalido
#define ARG_INV 105      // argumento lido da linha de comando invalido
#define FOPEN_ERR 106    // falha ao abrir arquivo
#define ZERO_DIV 107     // divisao por zero

// SIMD alignment macros
#define ALIGN_64 __attribute__((aligned(64)))
#define ALIGN_32 __attribute__((aligned(32)))
#define ALIGN_16 __attribute__((aligned(16)))

typedef double real_t;
typedef double rtime_t;
typedef char *string_t;

double timestamp(void);
void prnVetorFloat(float *x, int n);
void prnVetorDouble(double *x, int n);
void prnVetorLongDouble(long double *x, int n);

string_t markerName(string_t baseName, int n);
int isPot2(int n);

double *copyDoubleArray(double *a, int size); // copia array de double _a_ e retorna um ponteiro para copia
int max(int a, int b, int c);                 // retorna o maior de tres inteiros
char *getArgs(int argc, char **argv);         // separa os argumentos da linha de comando. retorna null caso nao haja arquivo de saida
void exitStatus(int code);                    // encerra o programa com status de saida diferente de 0
double sq_norma(double *array, int size);     // retorna o quadrado norma euclidiana do vetor _array_ de tamanho _size_
int isValidNum(double num);                   // testa se o numero nao e "nan" e nem "inf"
int pad(int n);                               // se n e potencia de 2, retorna n + 1; usado para evitar cache trashing

#endif

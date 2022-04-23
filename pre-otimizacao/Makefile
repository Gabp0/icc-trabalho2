# Gabriel de Oliveira Pontarolo GRR20203895
# Rodrigo Saviam Soffner GRR20205092

# Makefile

CC = gcc
CFLAGS = --std=c99 -Wall -O3 -mavx -march=native
CPPFLAGS = -I/usr/local/include 
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm -lmatheval
objects = main.o gaussianElimination.o luDecomposition.o gaussSeidel.o functions.o Rosenbrock.o utils.o newtonPadrao.o newtonModificado.o newtonInexato.o 

# default
all: newtonPC

# ligacao
newtonPC: $(objects)
	$(CC) -o newtonPC $(objects) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

# compilacao
main.o: main.c utils.h newtonPadrao.h newtonModificado.h
functions.o: functions.c functions.h utils.h
Rosenbrock.o: Rosenbrock.c Rosenbrock.h
gaussianElimination.o: gaussianElimination.c gaussianElimination.h utils.h 
luDecomposition.o: luDecomposition.c luDecomposition.h utils.h
gaussSeidel.o: gaussSeidel.c gaussSeidel.h utils.h
newtonInexato.o: newtonInexato.c newtonInexato.h Rosenbrock.h
newtonModificado.o: newtonModificado.c newtonModificado.h luDecomposition.h utils.h functions.h Rosenbrock.h
newtonPadrao.o: newtonPadrao.c newtonPadrao.h gaussianElimination.h utils.h functions.h Rosenbrock.h
utils.o: utils.c utils.h

# remove tudo que nao for codigo fonte
clean:
	-rm -f $(objects) *~ newtonPC
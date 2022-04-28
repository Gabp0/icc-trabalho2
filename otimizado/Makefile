# Gabriel de Oliveira Pontarolo GRR20203895
# Rodrigo Saviam Soffner GRR20205092

# Makefile

CC = gcc
CFLAGS = --std=c99 -Wall -O3 -mavx -march=native -fopt-info-vec-missed

LIKWID = /usr/bin/likwid
LIKWID_FLAGS = -I$(LIKWID)/include -DLIKWID_PERFMON
LIKWID_LIBS = -L$(LIKWID)/lib

CPPFLAGS = -I/usr/local/include $(LIKWID_FLAGS)
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm $(LIKWID_LIBS) -llikwid
objects = main.o gaussianElimination.o gaussSeidel.o functions.o Rosenbrock.o utils.o newtonPadrao.o newtonInexato.o 

# default
all: newtonPC

# ligacao
newtonPC: $(objects)
	$(CC) -o newtonPC $(objects) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

# compilacao
main.o: main.c utils.h newtonPadrao.h
functions.o: functions.c functions.h utils.h
Rosenbrock.o: Rosenbrock.c Rosenbrock.h
gaussianElimination.o: gaussianElimination.c gaussianElimination.h utils.h 
gaussSeidel.o: gaussSeidel.c gaussSeidel.h utils.h
newtonInexato.o: newtonInexato.c newtonInexato.h Rosenbrock.h
newtonPadrao.o: newtonPadrao.c newtonPadrao.h gaussianElimination.h utils.h functions.h Rosenbrock.h
utils.o: utils.c utils.h

# remove tudo que nao for codigo fonte
clean:
	-rm -f $(objects) *~ newtonPC

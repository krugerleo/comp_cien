CC=gcc
CFLAGS=-Wall -I/home/soft/likwid/include -DLIKWID_PERFMON -O3 -mavx -march=native
CLIBS=-lm -llikwid -L/usr/lib -L/home/soft/likwid/lib

cgSolver: cgSolver.c
	$(CC) $(CFLAGS) cgSolver.c -o cgSolver $(CLIBS)

clean:
	rm -rf *.txt *~ cgSolver

all: cgSolver

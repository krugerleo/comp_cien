CC=gcc
CFLAGS=-Wall -I/home/soft/likwid/include -DLIKWID_PERFMON -O0
CLIBS=-lm -llikwid -L/usr/lib -L/home/soft/likwid/lib

cgSolver: teste.c
	$(CC) $(CFLAGS) teste.c -o teste $(CLIBS)

clean:
	rm -rf *.txt *~ teste

all: teste

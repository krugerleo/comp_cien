CC=gcc
CFLAGS=-Wall -I/home/soft/likwid/include -DLIKWID_PERFMON -mavx -O3 -march=native
CLIBS=-lm -llikwid -L/usr/lib -L/home/soft/likwid/lib


teste: teste.c
	$(CC) $(CFLAGS) teste.c -o teste $(CLIBS)

clean:
	rm -rf *.txt *~ teste

all: teste

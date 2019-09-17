CC=gcc
CFLAGS=-Wall -I/home/soft/likwid/include -DLIKWID_PERFMON -mavx -O3 -march=native
CLIBS=-lm -llikwid -L/usr/lib -L/home/soft/likwid/lib


teste: main.c
	$(CC) $(CFLAGS) main.c -o main $(CLIBS)

clean:
	rm -rf *.txt *~ main

all: main

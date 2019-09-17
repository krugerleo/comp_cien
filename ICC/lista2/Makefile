CFLAGS = -O2 -std=c99 -lm

all: ex1 ex3 ex26

ex1: ex1.c 
	gcc -o ex1 ex1.c $(CFLAGS)

ex3: ex3.c
	gcc -o ex3 ex3.c $(CFLAGS)

ex26: ex26.c ../libmetodos/libZeroFuncao.c ../libmetodos/libPolinomios.c
	gcc -o ex26 ex26.c ../libmetodos/libZeroFuncao.c ../libmetodos/libPolinomios.c $(CFLAGS)

clean:
	rm -fr *.o

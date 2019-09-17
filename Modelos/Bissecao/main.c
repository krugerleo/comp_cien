#include <stdlib.h>
#include <stdio.h>

#define E 0.01

double bissecao( double(*f)(double), double a, double b ){
	double xn;
	do{
		xn = (a+b)/2;
		if(f(a)*f(xn) < 0)
			b = xn;
		else
			a = xn;
	}while( fabs(b-a) > E );
	
	return (b+a)/2;
}

double funcaoMatematica(double x){
	return (x*x - 4.0);
}

int main(){
	printf("Resultado da funcao matematica: %lf\n", bissecao(funcaoMatematica, 1.9,100000000000.0));
	return 0;
}

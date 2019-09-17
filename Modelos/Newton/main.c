#include <stdlib.h>
#include <stdio.h>

#define E 0.000000000000000001

double newton( double(*f)(double), double(*df)(double), double x0 ){
	double xk_velho;
	double xk_novo = x0;
	do{
		xk_velho = xk_novo;
		xk_novo = xk_velho - f(xk_velho)/df(xk_velho);
	}while( fabs(xk_novo - xk_velho) > E && fabs( f(xk_novo)-f(xk_velho) ) > E);
	return xk_novo;
}

double funcaoMatematica(double x){
	return (x*x*x - 2*x*x - 4*x + 3); //A raiz eh 0.1
}

double funcaoMatematica_derivada(double x){
	return (3*x*x - 4*x - 4);
}

int main(){
	printf("Resultado da funcao matematica: %lf\n", newton(funcaoMatematica,funcaoMatematica_derivada,-5));
	return 0;
}

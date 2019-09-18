#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

Polinomio pol_global;

double func_1 (double x)
{
  double px, dpx;
  calcPolinomioEDerivada(pol_global, x, &px, &dpx);
  return px;
}

double dfunc_1 (double x)
{
  double px, dpx;
  calcPolinomioEDerivada(pol_global, x, &px, &dpx);
  return dpx;
}

int main () {
	double px, dpx, raiz = 0.0;
	int it = 0;


	pol_global.grau = 2;
	pol_global.p[0] = -10.0;
	pol_global.p[1] = 0.0;
	pol_global.p[2] = 1;

	bisseccao (&func_1, 0, 5, 0.000000001, &it, &raiz);
	raiz = 0.0;
	it = 0;
	newton (&func_1, &dfunc_1, 5.0, 0.00000001, &it, &raiz);
	raiz = 0.0;
	it = 0;
	secante (&func_1, 0, 5, 0.00000001, &it, &raiz);

	return 0;
}


#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

#define __DEBUG__

int bisseccao (double (*f)(const double x), double a, double b,
               double eps, int *it, double *raiz)
{
  double erro, fa, fr, fb, raiz_ant;

#ifdef __DEBUG__
  fprintf(stderr, "\n\n****** Iterações Bissecção *******\n");
  fprintf(stderr, "=================================================================================================\n");
  fprintf(stderr, "  it            a            b           xs        f(xs)        |a-b|         f(a)         f(b)\n");
  fprintf(stderr, "=================================================================================================\n");
#endif

  	*it = 0;
  	do {
		raiz_ant = *raiz;
		*raiz = (a+b)/2;

		fa = f(a);
		fr = f(*raiz);

		if (fa*fr > 0)
			a = *raiz;
		else
			b  = *raiz;
		erro = fabs(*raiz - raiz_ant);
		*it = *it + 1;
	} while (erro > eps && *it < MAXIT);



#ifdef __DEBUG__
  fprintf(stderr, "%4d %12.8g %12.8g %12.8g ", *it, a, b, *raiz);
  fprintf(stderr, "%12.4e %12.4e %12.4e %12.4e\n", fr, fabs(a-b), fa, fb);
#endif

  return 0;
}


/**
 *
 */
int newton (double (*f)(const double x), double (*df)(const double x), double x0,
            double eps, int *it, double *raiz)
{
  double fx, dfx, erro;

#ifdef __DEBUG__
  fprintf(stderr, "\n\n****** Iterações Newton-Raphson *******\n");
  fprintf(stderr, "=================================================================================================\n");
  fprintf(stderr, "  it           x0           x0         raiz       f(raiz)   |raiz-x0|        f(x0)       df(x0)\n");
  fprintf(stderr, "=================================================================================================\n");
#endif

  	*it = 0;
	*raiz = x0;
  	do {
		x0 = *raiz;
		fx = f(x0);
		dfx = df(x0);
		*raiz = x0 - (fx/dfx);

		erro = fabs(*raiz - x0);
		*it = *it + 1;
	} while (erro > eps && *it < MAXIT);

#ifdef __DEBUG__
  fprintf(stderr, "%4d %12.8g %12.8g %12.8g ", *it, x0, x0, *raiz);
  fprintf(stderr, "%12.4e %12.4e %12.4e %12.4e\n", f(*raiz), fabs(*raiz-x0), fx, dfx);
#endif


  return 0;
}


/**
 *
 */
int secante (double (*f)(const double x), double x0, double x1,
             double eps, int *it, double *raiz)
{
  double fx0, fx1, erro, aux;
#ifdef __DEBUG__
  fprintf(stderr, "\n\n****** Iterações Secante *******\n");
  fprintf(stderr, "=================================================================================================\n");
  fprintf(stderr, "  it           x0           x1         raiz      f(raiz)    |raiz-x1|        f(x0)        f(x1)\n");
  fprintf(stderr, "=================================================================================================\n");
#endif

  	*it = 0;
  	do {
		aux = x1;
		fx0 = f(x0);
		fx1 = f(x1);

		x1 = x1 - ((fx1 * (x1 - x0)) / (fx1 - fx0));
		x0 = aux;

		erro = fabs(x0 - x1);
		*it = *it + 1;
	} while (erro > eps && *it < MAXIT);
	*raiz = x1;

#ifdef __DEBUG__
    fprintf(stderr, "%4d %12.8g %12.8g %12.8g ", *it, x0, x1, *raiz);
    fprintf(stderr, "%12.4e %12.4e %12.4e %12.4e\n", f(*raiz), fabs(*raiz-x0), fx0, fx1);
#endif

  return 0;
}

/**
 * Cálculo de Polinômios
 */
int calcPolinomioEDerivada(Polinomio pol, double x, double *px, double *dpx )
{
	double b;
	unsigned int i;

	*dpx = 0.0;
	b = pol.p[pol.grau];

	for (i = pol.grau - 1; (int)i >= 0; i--) {
		*dpx += b * pow (x, i);
		b = b * x + pol.p[i];
	}
	*px = b;

	return 0;
}

/**
 * Cálculo de Média
 */
double media(double *valores, unsigned long n)
{
	double soma = 0.0;
	unsigned long i = 0;

	for (i = 0; i < n; i++)
		soma += valores[i];

	return soma / n;

//	if (n == 2)
//	  return valores[0] + valores[1];
//	if (n == 1)
//	  return valores[0];
//
//	unsigned long div = n / 2;
//	return media(valores, div) + media(valores+div, n-div);
}



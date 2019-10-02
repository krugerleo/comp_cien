#include "utils.h"
#include <stdio.h>
#include <math.h>

// Retorna tempo em milisegundos
double timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

unsigned int encontraMaxPivo (SistLinear_t *SL, unsigned int diag) {
	unsigned int i, max = diag, n = SL->n;

	for (i = diag+1; i < n; i++)
		if (SL->A[i*n+diag] > SL->A[max*n+diag])
			max = i;
	return max;
}

void trocaLinha (SistLinear_t *SL, unsigned int linha1, unsigned int linha2) {
	unsigned int col, n = SL->n;
	real_t aux;
	for (col = 0; col < n; col++) {
		aux = SL->A[linha1*n+col];
		SL->A[linha1*n+col] = SL->A[linha2*n+col];
		SL->A[linha2*n+col] = aux;
		aux = SL->b[linha1];
		SL->b[linha1] = SL->b[linha2];
		SL->b[linha2] = aux;
	}
}

void retroSubs (SistLinear_t *SL, real_t *x) {
	real_t soma;
	int i = SL->n - 1, j, n = SL->n;

	x[i] = SL->b[i]/SL->A[i*n+i];
	for (i = n-2; i >= 0; i--) {
		soma = SL->b[i];
		for (j = i+1; j < n; ++j)
			soma -= SL->A[i*n+j] * x[j];
		x[i] = soma / SL->A[i*n+i];
	}
}

real_t calcNorma (real_t *x0, real_t *xk, int n) {
	int i;
	real_t norma = 0, aux;

	for (i = 0; i < n; i++) {
		aux = fabs (x0[i] - xk[i]);
		if (norma < aux)
			norma = aux;
	}
	return norma;
}

void zerarVetor (real_t *x, int n) {
	int i;
	for (i = 0; i < n; i++)
		x[i] = 0;

}

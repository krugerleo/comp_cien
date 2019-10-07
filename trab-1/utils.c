#include "utils.h"
#include <stdio.h>
#include <math.h>

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

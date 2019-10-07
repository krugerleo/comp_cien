#include "gaussSeidel.h"
#include "utils.h"
#include <stdio.h>

void gaussSeidel (sl5d *SL, double nl2r[], double tempos[], unsigned int maxIt) {
	unsigned int i, j, it, tam;

	tam = SL->nx * SL->ny;

	for (it = 0; it < maxIt; it++) {
		tempos[it] = timestamp();
		for (i = 0; i < tam; i++) {
			SL->x[i+3] = (SL->b[i] - (SL->ds[i]*SL->x[i+4] + SL->di[i]*SL->x[i+2] +
					SL->dsa*SL->x[i+6] + SL->dia*SL->x[i])) / SL->dp;
		}
		tempos[it] = timestamp() - tempos[it];
		nl2r[it] = calcNormaL2R(SL);
	}
}

double calcNormaL2R (sl5d *SL) {
	unsigned int i;
	double soma = 0.0;
	double r;

	for (i = 0; i < SL->tam; i++) {
		r = SL->b[i] - (SL->dp*SL->x[i+3] + SL->ds[i]*SL->x[i+4] + SL->di[i]*SL->x[i+2] +
				SL->dsa*SL->x[i+6] + SL->dia*SL->x[i]);
		r *= r;
		soma += r;
	}

	return sqrt(soma);
}

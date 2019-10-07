#include "gaussSeidel.h"
#include <stdio.h>

int main (int argc, char *argv[]) {
	sl5d SL;
	unsigned int i, j, nx, ny, maxIt;
	double *nl2r, *tempos;

	if (argc != 9) {
		printf ("erro!");
		return -1;
	}

	nx = atoi(argv[2]);
	ny = atoi(argv[4]);
	maxIt = atoi(argv[6]);

	nl2r = malloc(maxIt * sizeof(double));
	tempos = malloc(maxIt * sizeof(double));

	alocarSistLinear (&SL, nx, ny);

	gerarSistLinear (&SL);

	gaussSeidel (&SL, nl2r, tempos, maxIt);

	gerarSaidaGNUPLOT (&SL, nl2r, tempos, maxIt, argv[8]);

	desalocarSistLinear (&SL);
	free(nl2r);
	nl2r = NULL;
	free(tempos);
	tempos = NULL;

//	for (i = 0; i < SL.tam+6; i++)
//		printf("x%d = %.52f\n", i, SL.x[i]);

//	for (i = 0; i < maxIt; i++)
//		printf("r%d = %f\n", i, nl2r[i]);

//	for (i = 0; i < maxIt; i++)
//		printf("t%d = %5.3e\n", i, tempos[i]);

	return 0;
}

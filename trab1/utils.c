#include "utils.h"

// Retorna tempo em milisegundos
double timestamp(void) {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

int gerarSaidaGNUPLOT (sl5d *SL, double nl2r[], double tempos[], unsigned int tam, char *nomeArq) {
	FILE *arq;
	unsigned int i, j, k;
	double media, hx, hy, x, y;

	media = calcMedia (tempos, tam);
	arq = fopen(nomeArq, "w");

	if (arq == NULL)
		return -1;

	fprintf (arq,"###########\n# Tempo Método GS: %.15f\n#\n# Norma L2 do Residuo\n", media);

	for (i = 0; i < tam; i++) {

		fprintf (arq,"# i=%d: %.15f\n", i+1, nl2r[i]);
	}
	fprintf (arq,"###########\n#\n", media);
	fprintf (arq,"# Coluna 1 > x\n", media);
	fprintf (arq,"# Coluna 2 > y\n", media);
	fprintf (arq,"# Coluna 3 > u(x,y)\n#\n", media);

	hx = PI / (SL->nx + 1);
	hy = PI / (SL->ny + 1);

	k = 3;
	y = hy;

	for (j = 0; j < SL->ny; j++) {
		x = hx;
		for (i = 0; i < SL->nx; i++) {
			fprintf (arq,"%.15f %.15f %.15f\n", x, y, SL->x[k]);
			k++;
			x += hx;
		}
			y += hy;
	}

	fclose(arq);

	return 0;
}

double calcMedia (double v[], unsigned int tam) {
	double soma = 0.0, c = 0.0, y, t;
	unsigned int i;

	for (i = 0; i < tam; i++) {
		y = v[i] - c;
		t = soma + y;
		c = (t - soma) - y;
		soma = t;
	}
	return soma/tam;
}

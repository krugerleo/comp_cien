#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#define TAM 13
int main () {
	// inicializa gerador de nr aleatoreos
	srand(20192);

	SistLinear_t *SL;
	real_t *x = (real_t *) malloc(TAM * sizeof(real_t));

	SL = alocaSistLinear (TAM);
	inicializaSistLinear (SL, diagDominante, COEF_MAX);

	//SL = lerSistLinear();

	prnSistLinear (SL);

	//printf ("Iteracoes: %d\n",  gaussJacobi (SL, x, EPS));
	printf ("Iteracoes: %d\n",  gaussSeidel (SL, x, EPS));
	printf ("NormaL2 Residuo: %10g\n", normaL2Residuo(SL, x));

	prnVetor (x, TAM);

	(void) eliminacaoGauss(SL, x, 1);
	printf ("NormaL2 Residuo: %10g\n", normaL2Residuo(SL, x));
	prnVetor (x, TAM);


	liberaSistLinear (SL);
	free(x);
}


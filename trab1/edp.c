#include "edp.h"
#include <stdio.h>

/*!
	\brief Aloca memória para o sistema linear e calcula o tamanho da diagonal principal
	*
	\param SL Ponteiro para o sistema linear
	\param nx Número de pontos na dimensão x
	\param ny Número de pontos na dimensão y
*/
void alocarSistLinear (sl5d *SL, unsigned int nx, unsigned int ny) {
	SL->tam = nx*ny;
	SL->ds = malloc (SL->tam * sizeof(double));
	SL->di = malloc (SL->tam * sizeof(double));
	SL->b = malloc (SL->tam * sizeof(double));
	SL->x = malloc ((SL->tam+6) * sizeof(double));

	if (!(SL->ds) || !(SL->di) || !(SL->b) || !(SL->x)
			desalocarSistLinear(SL);

	SL->nx = nx;
	SL->ny = ny;
}

/*!
	\brief Desaloca memória do sistema linear
	*
	\param SL Ponteiro para o sistema linear
*/
void desalocarSistLinear (sl5d *SL) {
	free(SL->ds);
	free(SL->di);
	free(SL->b);
	free(SL->x);

}

/*!
	\brief Gera os coeficientes e  os termos independentes do sistema linear
	\brief Seta como 0 os valores do vetor SL->x[]
	*
	\param SL Ponteiro para o sistema linear
*/
void gerarSistLinear (sl5d *SL) {
	unsigned int i, j, k;
	double vDiagSup, vDiagInf, hx, hy, x, y;

	hx = PI / (SL->nx + 1);
	hy = PI / (SL->ny + 1);

	SL->dp = calcDiagPrin (hx, hy);
	vDiagSup = calcDiagSup (hx, hy);
	vDiagInf = calcDiagInf (hx, hy);
	SL->dsa = calcDiagSupAfas (hx, hy);
	SL->dia = calcDiagInfAfas (hx, hy);

	k = 0;
	y = hy;

	for (j = 0; j < SL->ny; j++) {
		x = hx;
		for (i = 0; i < SL->nx; i++) {
			SL->ds[k] = vDiagSup;
			SL->di[k] = vDiagInf;
			SL->b[k] = calcFxy(x, y);
			if (i == 0 && j > 0)
				SL->di[k] = 0;

			else if (i == (SL->nx-1) && j < (SL->ny-1))
				SL->ds[k] = 0;

			if (j == 0)
				SL->b[k] = SL->b[k] - calcUInf(x);
			else if (j == (SL->ny-1))
				SL->b[k] = SL->b[k] - calcUSup(x);
			k++;
			x += hx;
		}
			y += hy;
	}

	for (i = 0; i < SL->tam+6; i++)
		SL->x[i] = 0;
}

/*!
	\brief Calcula o valor da diagonal inferior
	*
	\param hx distância entre os pontos na dimensão x
	\param hy distância entre os pontos na dimensão y
*/
double calcDiagInf (double hx, double hy) {
	return (-2.0 - hx)/(2*(hx*hx));
//	return (hy*hy * (2.0*KAPA - LAMBDA*hx));
}

/*!
	\brief Calcula o valor da diagonal superior
	*
	\param hx distância entre os pontos na dimensão x
	\param hy distância entre os pontos na dimensão y
*/
double calcDiagSup (double hx, double hy) {
	return (-2.0 + hx)/(2*(hx*hx));
//	return (hy*hy * (2.0*KAPA + LAMBDA*hx));
}

/*!
	\brief Calcula o valor da diagonal inferior afastada
	*
	\param hx distância entre os pontos na dimensão x
	\param hy distância entre os pontos na dimensão y
*/
double calcDiagInfAfas (double hx, double hy) {
	return (-2.0 - hy)/(2*(hy*hy));
//	return (hx*hx * (2.0*KAPA - LAMBDA*hy));
}

/*!
	\brief Calcula o valor da diagonal superior afastada
	*
	\param hx distância entre os pontos na dimensão x
	\param hy distância entre os pontos na dimensão y
*/
double calcDiagSupAfas (double hx, double hy) {
	return (-2.0 + hy)/(2*(hy*hy));
//	return (hx*hx * (2.0*KAPA + LAMBDA*hy));
}

/*!
	\brief Calcula o valor da diagonal principal
	 *
	\param hx distância entre os pontos na dimensão x
	\param hy distância entre os pontos na dimensão y
*/
double calcDiagPrin (double hx, double hy) {
	return (2.0 * (hx*hx + hy*hy) + (4.0*PI*PI)*(hx*hx*hy*hy))/(hx*hx*hy*hy);
//	return (-4.0 * (hx*hx + hy*hy) + ETA*2.0*hx*hx*hy*hy);
}


/*!
	\brief Calcula o valor da função f(x,y) no ponto x e y passados por paramêtro.
	 *
	\param x valor do ponto x
	\param y valor do ponto y
*/
double calcFxy (double x, double y) {
	return 4.0*PI*PI * (sin(2.0*PI*x) * sinh(PI*y) + sin(2.0*PI*(PI-x)) * sinh(PI*(PI-y)));
}

/*!
	\brief Calcula o valor da função u(x,0) no ponto x passado por paramêtro.
	 *
	\param x valor do ponto x
*/
double calcUInf (double x) {
	return (sin(2.0*PI*(PI-x)) * sinh(PI*PI));
}

/*!
	\brief Calcula o valor da função u(x,pi) no ponto x passado por paramêtro.
	 *
	\param x valor do ponto x
*/
double calcUSup (double x) {
	return (sin(2.0*PI*x) * sinh(PI*PI));
}

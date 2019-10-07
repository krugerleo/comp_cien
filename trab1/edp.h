#ifndef __EDP__
#define __EDP__

#include <math.h>
#include <stdlib.h>
#define PI 3.14159265358979323846
#define KAPA -1.0
#define LAMBDA 1.0
#define ETA 4*PI*PI

typedef struct t_SistLinear5Diag {
	double dp, *ds, dsa, *di, dia, *b, *x;
	unsigned int nx, ny, tam;
} sl5d;

void alocarSistLinear (sl5d *SL, unsigned int nx, unsigned int ny);

void desalocarSistLinear (sl5d *SL);

void gerarSistLinear (sl5d *SL);

double calcDiagInf (double hx, double hy);

double calcDiagSup (double hx, double hy);

double calcDiagInfAfas (double hx, double hy);

double calcDiagSupAfas (double hx, double hy);

double calcDiagPrin (double hx, double hy);

double calcFxy (double x, double y);

double calcUInf (double x);

double calcUSup (double x);

#endif


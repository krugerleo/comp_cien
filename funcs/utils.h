#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <sys/time.h>
#include "SistemasLineares.h"

double timestamp(void);
unsigned int encontraMaxPivo (SistLinear_t *SL, unsigned int diag);
void trocaLinha (SistLinear_t *SL, unsigned int linha1, unsigned int linha2);
void retroSubs (SistLinear_t *A, real_t *x);
real_t calcNorma (real_t *x0, real_t *xk, int n);
void zerarVetor (real_t *x, int n);

#endif // __UTILS_H__


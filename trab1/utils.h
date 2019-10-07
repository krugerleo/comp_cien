#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "edp.h"

double timestamp(void);

int gerarSaidaGNUPLOT (sl5d *SL, double nl2r[], double tempos[], unsigned int tam, char *nomeArq);

double calcMedia (double v[], unsigned int tam);

#endif // __UTILS_H__


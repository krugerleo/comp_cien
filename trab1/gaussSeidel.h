#ifndef __GAUSSSEIDEL__
#define __GAUSSSEIDEL__

#include "edp.h"
#include "utils.h"

void gaussSeidel (sl5d *SL, double nl2r[], double tempos[], unsigned int maxIt);

double calcNormaL2R (sl5d *SL);

#endif

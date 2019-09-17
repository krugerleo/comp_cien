#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libPolinomios.h"


int bissecao(double *p, int n, double *x, double a, double b, double erroMax) {
    double fa, fb, fx, erro, xVelho;
    fa = px(p, n, a);
    fb = px(p, n, b);

    do {
        xVelho = *x;
        *x = (a+b)/2;
        fx = px(p, n, *x);
        if((fa * fx) < 0) {
            b = *x;
        } else {
            a = *x;
        }
        erro = fabs(*x - xVelho);
    } while((fa * fb) != 0.0f && erro > erroMax);
    return 1;
}

int newtonRaphson(double *p, int n, double *x, double erroMax) {
    double px, dpx, erro, xNovo;

    do {
        pxDpx(p, n, *x, &px, &dpx);
        if(dpx == 0.0f) return -1;
        xNovo = *x - px/dpx;
        erro = fabs(xNovo - *x);
        *x = xNovo;
    } while(erro > erroMax);
    return 0;
}

int secante(double *p, int n, double *x, double x0, double erroMax) {
    double pk, pxk, xNovo, erro;

    do {
        pk = px(p, n, *x);
        pxk = px(p, n, x0);
        xNovo = *x - (pk * (*x - x0)) / (pk - pxk);
        erro = fabs(xNovo - *x);
        *x = xNovo;
    } while(erro > erroMax);
    return 0;
}

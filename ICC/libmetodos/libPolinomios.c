#include <stdio.h>
#include <stdlib.h>

double px(double *p, int n, double x) {
    double px = p[n];

    for(int i = n-1; i >= 0; --i) {
        px = p[i] + px*x;
    }
    return px;
}

void pxDpx(double *p, int n, double x, double *px, double *dpx) {
    double b = p[n], c = p[n];

    for(int i = n-1; i > 0; --i) {
        b = p[i] + b*x;
        c = b + c*x;
    }
    b = p[0] + b*x;
    *px = b;
    *dpx = c;
}

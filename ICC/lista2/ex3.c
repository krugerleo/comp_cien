/* Calcular a raiz de f(x) = x² + x - 6 com erro relativo < 0.1
 * a) pela bisseção com xl = 0, xu = 3
 * b) por newton-raphson com x0 = -1.5 e depois x0 = 1.5
 * c) pela secante utilizando x0 = 1.5, x1 =1.8
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * Calcula f(x)
 */
double px(double *p, int n, double x) {
    double px = p[n];

    for(int i = n-1; i >= 0; --i) {
        px = p[i] + px*x;
    }
    return px;
}

void calculaPolinomioEDerivada(double *p, int n, double x, double *px, double *dpx) {
    double b = p[n], c = p[n];
    for(int i = n-1; i > 0; --i) {
        b = p[i] + b * x;
        c = b + c * x;
    }
    b = p[0] + b * x;
    *px = b;
    *dpx = c;
}

int bissecao(double *p, int n, double *x, double a, double b, double erroMax) {
    double fa, fb, fx, erro, xOld;
    fa = px(p, n, a);
    fb = px(p, n, b);
    
    do {
        xOld = *x;
        *x = (a+b)/2;
        fx = px(p, n, *x);
        if(fa*fx < 0) {
            b = *x;
        } else {
            a = *x;
        }
        // Erro relativo
        erro = fabs(*x - xOld);
    }while(fa*fb != 0.0 && erro > erroMax);
    return 1;
}

int newtonRaphson(double *p, int n, double *x, double erroMax) {
    double px, dpx, erro, xNovo;

    do {
        calculaPolinomioEDerivada(p, n, *x, &px, &dpx);
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

int main() {
    double *p, x, erro = 0.001;
    int n;

    printf("Grau do polinomio: \n");
    scanf("%d", &n);

    p = malloc((n+1) * sizeof(double));

    for(int i = 0; i <= n; ++i){
        printf("x^%d\t", i);
        scanf("%lf", &p[i]);
    }
    puts("");
    
    x = 0;
    bissecao(p, n, &x, 0.0, 3.0, erro);
    printf("Bissecao: x = %.20lf\n", x);

    x = -1.5f;
    newtonRaphson(p, n, &x, erro);
    printf("Newton (-1.5): x = %.20lf\n", x);
    x = 1.5f;
    newtonRaphson(p, n, &x, erro);
    printf("Newton (1.5): x = %.20lf\n", x);

    x = 1.8f;
    secante(p, n, &x, 1.5f, erro);
    printf("Secante: x = %.20lf\n", x);

    return 0;
}

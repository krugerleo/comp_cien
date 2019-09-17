#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Essa função é na verdade o método de Newton-Raphson cuja função de iteração
// é i(x) = x - f(x)/f'(x);
/**
 *  p: coeficientes de um polinomio
 *  n: grau do polinomio p
 *  x: valor inicial e resposta
 *  erroMax: maior erro aceitavel
 */
int funcaoFazAlgo(double *p, int n, double *x, double erroMax) {
    double px, dpx, erro, x_new;
    do {
        calculaPolinomioEDerivada(p, n, *x, &px, &dpx);
        // Comparar usando == não é uma comparação segura em caso de ponto flutuante
        if(dpx == 0.0f) return -1;
        x_new = *x - px / dpx;
        erro = fabs(x_new - *x);
        *x = n_new;
    } while(erro > erroMax);
    return 0;
}

/**
 *  p, n, x: idem a funcão anterior
 *  px: valor do polinomio no ponto x
 *  dpx: valor da primeira derivada no ponto x
 */
void calculaPolinomioEDerivada(double *p, int n, double x, double *px, double *dpx) {
    
}

int main() {
    return 0;
}

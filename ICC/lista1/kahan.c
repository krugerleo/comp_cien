#include <stdio.h>
#include <stdlib.h>

double somaKahanDouble(double *dados, unsigned int tam) {
    double soma = 0.0;
    double c = 0.0;
    for(unsigned int i=0; i < tam; ++i) {
        double y = dados[i] - c;
        double t = soma + y;
        c = (t - soma) - y;
        soma = t;
    }
    return soma;
}

float somaKahanFloat(float *dados, unsigned int tam) {
    float soma = 0.0;
    float c = 0.0;
    for(unsigned int i=0; i < tam; ++i) {
        float y = dados[i] - c;
        float t = soma + y;
        c = (t - soma) - y;
        soma = t;
    }
    return soma;
}

#include <stdio.h>
#include <stdlib.h>
#include "kahan.h"

#define VALOR 0.6f
#define NUM_ELEMENTOS 1000000

float somaSequencial(float *dados, unsigned int tam) {
    float soma = 0.0;
    
    while(tam--) {
        soma += dados[tam];
    }
    // O return não está no código da lista, mas faz sentido tê-lo
    return soma;
}

float somaPar(float *dados, unsigned int tam) {
    if(tam == 2) return dados[0] + dados[1];
    if(tam == 1) return dados[0];
    
    unsigned int div = tam / 2;
    return somaPar(dados, div) + somaPar(dados, tam-div);
}

int main() { 
    // Preenche um vetor
    float *dados = malloc(NUM_ELEMENTOS * sizeof(float));
    for(unsigned int i=0; i < NUM_ELEMENTOS; i++) dados[i] = VALOR;

    float soma1 = somaSequencial(dados, NUM_ELEMENTOS);
    printf("Soma sequencia: %1.15f\n", soma1);

    float soma2 = somaPar(dados, NUM_ELEMENTOS);
    printf("Soma par: %1.15f\n", soma2);

    float soma3 = somaKahanFloat(dados, NUM_ELEMENTOS);
    printf("Soma Kahan: %1.15f\n", soma3);

    free(dados);

    return 0;
}

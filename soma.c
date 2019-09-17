#include <stdio.h>
#include <stdlib.h>
#define VALOR         0.6f
#define NUM_ELEMENTOS 1000
#include <stdint.h>
#include <float.h>
#include <math.h>
typedef union
{
    int32_t i;
    float f;
    struct
    {   // Bitfields for exploration. Do not use in production code.
        uint32_t mantissa : 23;
        uint32_t exponent : 8;
        uint32_t sign : 1;
    } parts;
} Float_t;
void printFloat_t( Float_t num )
{
   printf("f:%1.9e, ix:0x%08X, s:%d, e:%d, mx:0x%06X\n",
            num.f, num.i,
            num.parts.sign, num.parts.exponent, num.parts.mantissa);
}
float somaSequencia( float *dados, unsigned int tam )
{
    float soma = 0.0;
    while ( tam-- )
    {
        soma += dados[tam];
    }
}
float somaPar( float *dados, unsigned int tam )
{
    if (tam == 2)
        return dados[0] + dados[1];
    if (tam == 1)
        return dados[0];
    unsigned int div = tam / 2;
    return somaPar(dados, div) + somaPar(dados+div, tam-div);
}
float somaKahan( float *dados, unsigned int tam )
{
    float soma = 0.0;
    float c = 0.0;
    float x,y;
    while( tam-- ){
        y = dados[tam] - c;
        x = soma + y;
        c = (x - soma) - y;
        soma = x;
    }
    return soma;
}
void main()
{
    Float_t x;
    // Preenche um vetor
    float *dados = (float*) malloc(NUM_ELEMENTOS * sizeof(float));
    for (unsigned int i = 0; i < NUM_ELEMENTOS; ++i)
        dados[i] = VALOR;
    float soma1 = somaSequencia( dados, NUM_ELEMENTOS );
    printf("Soma sequencia: %1.15f\n", soma1);
    float soma2 = somaPar( dados, NUM_ELEMENTOS );
    printf("Soma par: %1.15f\n", soma2);
    float soma3 = somaKahan( dados, NUM_ELEMENTOS );
    printf("Soma kahan: %1.15f\n", soma3);
    x.i = 47;
    printFloat_t(x);
    free (dados);
}

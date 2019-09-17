#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

typedef union {
    int32_t i;
    float f;
    // Bitfields
    struct {
        uint32_t mantissa: 23;
        uint32_t exponent: 8;
        uint32_t sign: 1;
    } parts;
}Float_t;

void exploreFloat(Float_t ft) {
    printf("Float: %1.15e\t Int: %i\t Mantissa: %i\t Expoente: %i\t Sinal:%i\n",
            ft.f, ft.i, ft.parts.mantissa, ft.parts.exponent, ft.parts.sign);
}

void nextFive(float f) {
    Float_t ft;
    ft.f = f;
    int i = 5;
    while(i--) {
        f = ft.f;
        ft.parts.mantissa += 1;
        exploreFloat(ft);
        float difference = fabsf(ft.f - f);
        printf("Diferença entre %1.15e e %1.15e é %1.15e\n", ft.f, f, difference);
    }
}

int main() {
    float f;
    scanf("%f", &f);
    nextFive(f);
    return 0;
}

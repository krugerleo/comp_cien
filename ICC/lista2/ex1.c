#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// Função original
bool baskara(double a, double b, double c, double *r) {
    double delta, ml;
    delta = b*b - 4*a*c;
    // Comparação com 0.0f não funciona sempre
    // Por causa de possíveis erros causados na subtração, o delta pode dar
    // um valor muito próximo de 0 quando deveria dar 0
    if(delta < 0.0f) return false;

    ml = sqrt(delta);
    r[0] = (-b + ml) / (2.0  * a);
    r[1] = (-b - ml) / (2.0  * a);
    return true;
}

int main() {
    double r[2];
    double a, b, c;
    scanf("%lf %lf %lf", &a, &b, &c);
    if(baskara(a, b, c, r)) 
        printf("%1.15e %1.15e\n", r[0], r[1]);
    else
        printf("Nao ha raizes\n");
    return 0;
}

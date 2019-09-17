#include <stdio.h>
#include <stdlib.h>
#include "../libmetodos/libZeroFuncao.h"

int main() {
    int n = 3;
    double p[4] = {3.993E-4, 0, -0.165, 1}, x;
    
    x = 1;
    newtonRaphson(p, n, &x, 0.0001);
    printf("Zero de função de x³ - 0.165x² + 3.993*10⁻⁴: %.20lf\n", x);
    return 0;
}

/**
 * Calcula o valor de x dado um polinômio p de grau n
 */
double px(double *p, int n, double x);

/**
 * Calcula o valor de x dado um polinômio p de grau n e sua derivada no ponto
 */
void pxDpx(double *p, int n, double x, double *px, double *dpx);

/**
 * Método da bisseção
 * Calcula o zero de função de um polinômio p de grau n buscando no intervalo (a,b)
 * parando quando o erro for menor do que erroMax
 */
int bissecao(double *p, int n, double *x, double a, double b, double erroMax);

/**
 * Método de Newton-Raphson
 * Calcula o zero de função usando como função de iteração i(x) = x - f(x)/f'(x)
 */
int newtonRaphson(double *p, int n, double *x, double erroMax);

/**
 * Método da Secante
 * Calcula o zero de função trocando a derivada do método de Newton por diferenças finitas
 */
int secante(double *p, int n, double *x, double x0, double erroMax);

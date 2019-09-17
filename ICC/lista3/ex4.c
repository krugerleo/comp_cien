#include <stdio.h>
#include <stdlib.h>

void eliminaGauss(double **a, int n, double *b) {
    for(int k=0; k < n; ++k) {
        for(int i = k+1; i < n; ++i) {
            double m = a[i][k]/a[k][k];
            a[i][k] = 0.0f;
            for(int j=k+1; j < n; ++j) {
                a[i][j] -= a[k][j] * m;
            }
            b[i] -= b[k] * m;
        }
    }
}

void eliminaGaussPivo(double **a, int n, double *b) {
    for(int k=0; k < n; ++k) {
        int max = k;
        for(int i = k+1; i < n; ++i) {
            max = (a[i][k] > a[max][k]) ? i : max;
        }
        if(max != k) trocaLinha(a, n, b, max, k);
        for(int i = k+1; i < n; ++i) {
            double m = a[i][k]/a[k][k];
            a[i][k] = 0.0f;
            for(int j=k+1; j < n; ++j) {
                a[i][j] -= a[k][j] * m;
            }
            b[i] -= b[k] * m;
        }
    }
}

void trocaLinha(double **a, int n, double *b, int l1, int l2) {
    double temp;
    temp = b[l2];
    b[l2] = b[l1];
    b[l1] = temp;
    for(int i = 0; i < n; ++i) {
        temp = a[l1][i];
        a[l1][i] = a[l2][i];
        a[l2][i] = temp;
    }
}

void retro(double **a, int n, double *b, double *x) {
    x[n-1] = b[n-1]/a[n-1][n-1];
    for(int i = n-2; i >= 0; --i) {
        double soma=0.0f;
        for(int j = n-1; j > i; --j) {
            soma += a[i][j]*x[j];
        }
        x[i] = (b[i] - soma)/a[i][i];
    }
    
}

int main() {
    double **a, **A, *b, *B, *x;
    int n;

    printf("Ordem da matriz: \n");
    scanf("%d", &n);

    a = malloc(n*sizeof(double*));
    A = malloc(n*sizeof(double*));
    b = malloc(n*sizeof(double));
    B = malloc(n*sizeof(double));
    x = malloc(n*sizeof(double));
    for(int i=0; i < n; ++i) {
        a[i] = malloc(n*sizeof(double));
        A[i] = malloc(n*sizeof(double));
    }

    for(int i=0; i < n; ++i) {
        for(int j=0; j < n; ++j) {
            printf("(%d, %d)\t", i, j);
            scanf("%lf", &a[i][j]);
            A[i][j] = a[i][j];
        }
    }

    printf("B\n");
    for(int i=0; i < n; ++i) {
        printf("b[%d]", i);
        scanf("%lf", &b[i]);
        B[i] = b[i];
    }

    eliminaGauss(a, n, b);
    retro(a, n, b, x);
    printf("Resposta: \n");
    for(int i = 0; i < n; ++i) {
        printf("%.4lf\n", x[i]);
    }

    eliminaGaussPivo(A, n, B);
    retro(A, n, B, x);
    printf("Pivo: \n");
    for(int i = 0; i < n; ++i) {
        printf("%.4lf\n", x[i]);
    }

    return 0;
}


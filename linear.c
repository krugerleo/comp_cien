#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N 3

void imprimeMatriz(double** A){
	int i, j;
	for(i = 0; i < N; i++){
		printf("|");
		for(j = 0; j < N; j++){
			printf(" %f", A[i][j]);
		}
		printf(" |\n");
	}
}

void imprimeVetor(double* v){
	int i;
	for(i = 0; i < N; i++){
		printf("|%f|\n",v[i]);
	}
}
void retrosub(double **dados,double *valores,unsigned int tam){
    double *incog = malloc(sizeof(tam)); // vetor valores de incog
    for(int x = tam-1 ; x >= 0; x--){
        incog[x] = valores[x];
        for(unsigned int j = tam-1; j > x; j--){
            incog[x] -= dados[x][j] * incog[j];
        }
        incog[x] /= dados[x][x];
        printf(" = %+.2f\n", incog[x]);
    }
}

double** geraMatrizEscalonada(double* x, double* b){

	double** A = malloc(sizeof(double*)*N);
	srand(time(NULL));

	for(unsigned int i = 0; i < N; i++){
		A[i] = malloc(sizeof(double)*N);
	}

	A[0][0] = 1; A[0][1] = 2; A[0][2] = 3; b[0] = 10;
	A[1][0] =  0; A[1][1] = 2; A[1][2] = 3; b[1] = 7;
	A[2][0] = 0; A[2][1] =  0; A[2][2] = 3; b[2] = 3;

	printf("Vetor b: \n");
	imprimeVetor(b);
	printf("Matiz A: \n");
	imprimeMatriz(A);
	return A;
}
void gauss(float dados[3][3],float *valores,unsigned int tam){

}

int main() {
    double *x,*y;
    double **mtx = geraMatrizEscalonada(x,y);
    retrosub(mtx,y,N);
}
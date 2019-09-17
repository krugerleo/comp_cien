#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define N 10

double staticX[N] = {10000.2 ,20000.3 ,30000.4, 40000.5, 50000.6 ,60000.7, 70000.8, 80000.9, 90000.10, 100000.11};

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
double** geraMatrizEscalonada(double* x, double* b){
	double** A = malloc(sizeof(double*)*N);
	srand(time(NULL));
	
	int i, j;
	
	for(i = 0; i < N; i++){
		b[i] = 0.0;
		A[i] = malloc(sizeof(double)*N);
		for(j = 0; j < i; j++){
			A[i][j] = 0.0;
		}
		
		while(j < N){
			A[i][j] = 10*(double)rand()/(double)RAND_MAX;
			b[i]   += A[i][j]*x[j];
			j++;
		}
	}
	
	printf("Vetor b: \n");
	imprimeVetor(b);
	printf("Matiz A: \n");
	imprimeMatriz(A);
	return A;
}

void retroSubstituicao(double** A, double* x, double* b, int n){
	int i, j;
	x[n-1] = b[n-1]/A[n-1][n-1];
	for( i = n-2; i >=0; --i){
		x[i] = b[i];
		for(j = n-1; j > i; --j){
			x[i] -= A[i][j]*x[j];
		}
		x[i] /= A[i][i];
	}
}



int main(){
	double* b = malloc(sizeof(double)*N);
	double** A = geraMatrizEscalonada(staticX,b);
	
	double* x = malloc(sizeof(double)*N);
	retroSubstituicao(A,x,b,N);
	
	printf("\n Resultado: \n");
	imprimeVetor(x);
	
	return 0;
}

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define N 4

double staticX[N] = {10000.2 ,20000.3 ,30000.4, 40000.5};

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
		A[i] = malloc(sizeof(double)*N);
	}


	A[0][0] = -1; A[0][1] = -4; A[0][2] = 2; A[0][3] =  1; b[0] = -32;
	A[1][0] =  2; A[1][1] = -1; A[1][2] = 7; A[1][3] =  9; b[1] = 14;
	A[2][0] = -1; A[2][1] =  1; A[2][2] = 3; A[2][3] =  1; b[2] = 11;
	A[3][0] =  1; A[3][1] = -2; A[3][2] = 1; A[3][3] = -4; b[3] = -4;

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

void trocaLinhas(double** A, double*b, int k, int l){
	int i;
	double aux;
	for(i = 0; i < N; i++){
		aux = A[k][i];
		A[k][i] = A[l][i];
		A[l][i] = aux;
	}
	aux = b[k];
	b[k] = b[l];
	b[l] = aux;
}

void metodoDeGauss(double** A, double*b, int n){
	int j, k, i, l;
	double m;

	for(j=0; j<n-1; j++){

		//Pivotamento
		k = j;
		for(i = j+1; i < n; i++){
			if( abs(A[i][j]) > abs(A[k][j])){
				k = i;
			}
		}
		trocaLinhas(A,b,k,j);
		//
		imprimeMatriz(A);
		imprimeVetor(b);
		puts("");

		for(i = j+1; i < n; i++){

			m = A[i][j]/A[j][j];
			A[i][j] = 0;
			printf(">%f\n",m);
			for(l = j+1; l < n; l++){
				A[i][l] = A[i][l] - m*A[j][l];
			}
			b[i] = b[i] - m*b[j];
		}
		//
		imprimeMatriz(A);
		imprimeVetor(b);
		puts("");
	}
}



int main(){
	double* b = malloc(sizeof(double)*N);
	double** A = geraMatrizEscalonada(staticX,b);

	double* x = malloc(sizeof(double)*N);

	puts("----------GAUS-----------");
	metodoDeGauss(A,b,N);
	puts("-------------------------");
	puts("Apos Gauss:");
	imprimeMatriz(A);
	imprimeVetor(b);

	retroSubstituicao(A,x,b,N);

	printf("\n Resultado: \n");
	imprimeVetor(x);

	return 0;
}

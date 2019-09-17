#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define N 4

double staticX[N] = {10000.2 ,20000.3 ,30000.4, 40000.5};

void imprimeMatriz(double** A, int n){
	int i, j;
	for(i = 0; i < n; i++){
		printf("|");
		for(j = 0; j < n; j++){
			printf(" %.3f", A[i][j]);
		}
		printf(" |\n");
	}
}

void imprimeVetor(double* v,int n){
	int i;
	for(i = 0; i < n; i++){
		printf("|%f|\n",v[i]);
	}
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

void trocaLinhas(double** A, double*b, int k, int l, int n){
	int i;
	double aux;
	for(i = 0; i < n; i++){
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
		trocaLinhas(A,b,k,j,n);
		//
		
		for(i = j+1; i < n; i++){
			
			m = A[i][j]/A[j][j];
			A[i][j] = 0;
			for(l = j+1; l < n; l++){
				A[i][l] = A[i][l] - m*A[j][l];
			}
			b[i] = b[i] - m*b[j];
		}
		//
	}
}

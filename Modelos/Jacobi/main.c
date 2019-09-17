#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 3
#define maxIte 10000000
#define erroMax 0.00001

double staticX[N] = {10000.2 ,20000.3 ,30000.4};

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
		/*
		b[i] = 0.0;
		
		for(j = 0; j < i; j++){
			A[i][j] = 0.0;
		}
		
		while(j < N){
			A[i][j] = 10*(double)rand()/(double)RAND_MAX;
			b[i]   += A[i][j]*x[j];
			j++;
		}
		*/
	}
	
	A[0][0] = 12; A[0][1] = 3; A[0][2] = -5; b[0] =  1;
	A[1][0] =  1; A[1][1] = 5; A[1][2] =  3; b[1] = 28;
	A[2][0] =  3; A[2][1] = 7; A[2][2] = 13; b[2] = 76;
	
	/*
	A[0][0] = -1; A[0][1] = -4; A[0][2] = 2; A[0][3] =  1; b[0] =  -32;
	A[1][0] =  2; A[1][1] = -1; A[1][2] = 7; A[1][3] =  9; b[1] =  14;
	A[2][0] = -1; A[2][1] =  1; A[2][2] = 3; A[2][3] =  1; b[2] =  11;
	A[3][0] =  1; A[3][1] = -2; A[3][2] = 1; A[3][3] = -4; b[3] =  -4;
	*/
	printf("Vetor b: \n");
	imprimeVetor(b);
	printf("Matiz A: \n");
	imprimeMatriz(A);
	return A;
}

void forwardSubstitution(double** A, double* x, double* b, int n){
	int i,j;
	x[0] = b[0] / A[0][0];
	for(i = 1; i < n; i++){
		x[i] = b[i];
		for(j = 0; j < i; j++){
			x[i] -= A[i][j]*x[j];
			
		}
		x[i] /= A[i][i];
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

double maxDif(double** x, int n){
	double maior = fabs(x[0][0] - x[1][0]);
	double aux;
	int i;
	for(i = 1; i < n; i++){
		aux = fabs(x[0][i] - x[1][i]);
		if(aux > maior)
			maior = aux;
	}
	return maior;
}


double* jacobi(double** A, double* b, double* x1,int n){
	double* x2 = malloc(sizeof(double)*n);
	double** x = malloc(sizeof(double*)*2);
	x[0] = x1;
	x[1] = x2;
	int i,j;
	
	int ite = 0;
	int kAtual = 0;
	int kAntigo = 1;
	
	for(i=0;i<n;i++)
		x[0][i] = 0.0;
	
	double somatorio;
	do{
		if(kAtual == 0){
			kAtual = 1;
			kAntigo = 0;
		}
		else{
			kAtual = 0;
			kAntigo = 1;
		}
	
		for(i=0; i < n; i++){
			somatorio = 0.0;
			for(j = 0; j < n; j++){
				if(j != i){
					somatorio += (A[i][j]*x[kAntigo][j])/A[i][i];
				}
			}
			x[kAtual][i] = (b[i]/ A[i][i] - somatorio) ;
		}
		/*
		printf("--Atual-- %d\n",kAtual);
		imprimeVetor(x[kAtual]);
		printf("---------\n--Antigo-- %d\n",kAntigo);
		imprimeVetor(x[kAntigo]);
		puts("---------\n");
		*/
		ite++;
	} while(ite<=maxIte && maxDif(x,n) > erroMax);
	printf(">Resolveu em %d iteracoes!\n",ite);
	
	return x[kAtual];
}


int main(){
	double* b = malloc(sizeof(double)*N);
	double** A = geraMatrizEscalonada(staticX,b);
	
	
	double* x = malloc(sizeof(double)*N);
	
	
	
	puts("----------Jacobi-----------");
	double* y = jacobi(A,b,x,N);
	puts("-------------------------");
	puts("Apos Jacobi:");
	puts("x:");
	imprimeVetor(y);
	return 0;
}

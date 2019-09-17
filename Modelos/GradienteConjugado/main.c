#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 3
#define maxIte 10000000
#define erroMax 0.00001

//double staticX[N] = {10000.2 ,20000.3 ,30000.4};


/***********************
 * N: tamanho do sistema linear
 * k: numero da diagonal, 0 = diagonal principal, 1 = acima/abaixo da diagonal, 2 = ...
 * kMax: numero de bandas do sistema linear
 * diag: vetor para armazenar os valores da diagonal. Deve ser alocado por quem chama a função.
 ***********************/
int generateRandomDiagonal( unsigned int n, unsigned int k, unsigned int kMax, double *diag )
{
  if ( !diag || n < 3 || kMax > n/2 || k < 0 || k > kMax )
    return (-1);

  /* garante valor dominante para diagonal principal */
  double fator;// = (k == 0) ? ((double)(kMax-1)) : (0.0);
  if(k == 0)
	fator = (double)(kMax);
  else
	fator = 0.0;
  double invRandMax = 1.0 / (double)RAND_MAX;
  int i;
  for (i=0; i < (n - k) ; ++i)
  {
    diag[i] = fator + (double)rand() * invRandMax;
  }

  return (0);
}




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
double** geraMatrizEscalonada(double* b){
	
	double** A = malloc(sizeof(double*)*N);
	srand(time(NULL));
	
	
	int i, j;
	
	for(i = 0; i < N; i++){
		A[i] = malloc(sizeof(double)*N);
		for(j = 0; j < N; j++)
			A[i][j] = 0.0;
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
	
	
	//f(x) = 4π² ( sin(2πx) + sin(2π(π-x)) )
	double f(double xi){
		double auxFloat1 = 2*M_PI*xi;
		double auxFloat2 = 2*M_PI*(M_PI-xi);
		
		return 4*M_PI*M_PI * ( sin(auxFloat1) + sin(auxFloat2) );
	}
	
	double* aux = malloc(sizeof(double));
	generateRandomDiagonal(N,0,1,aux);
	for(i = 0; i < N; i++){
		A[i][i] = aux[i];
		b[i] = f( (double) (i * M_PI / N) );// para todo i = 0:n-1, 
	}
	generateRandomDiagonal(N,1,1,aux);
	for(i = 0; i < N; i++){
		if(i+1 >= 0 && i+1 < N)
			A[i][i+1] = aux[i];
	}
	generateRandomDiagonal(N,1,1,aux);
	for(i = 0; i < N; i++){
		if(i-1 >= 0 && i-1 < N)
			A[i][i-1] = aux[i];
	}
	
	/*
	A[0][0] = 12; A[0][1] = 3; A[0][2] = -5; b[0] =  1;
	A[1][0] =  1; A[1][1] = 5; A[1][2] =  3; b[1] = 28;
	A[2][0] =  3; A[2][1] = 7; A[2][2] = 13; b[2] = 76;
	*/
	
	/*
	A[0][0] = 2; A[0][1] = 1; A[0][2] = 2; b[0] =  4;
	A[1][0] = 1; A[1][1] = 2; A[1][2] = 1; b[1] = -1;
	A[2][0] = 3; A[2][1] = 5; A[2][2] = 2; b[2] =  1;
	*/
	
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

void multiplicaMatrizPorVetor(double** A, double* v, double* r, int n){
	int i,j;
	for(i = 0; i < n; i++){
		r[i] = 0.0;
		for(j = 0; j < n; j++){
			r[i] += A[i][j]*v[j];
		}
	}
}


double multiplicaVetores(double* v1, double* v2, int n){
	double soma = 0.0;
	int i;
	for(i = 0; i < n; i++){
		soma += (v1[i]*v2[i]);
	}
	return soma;
}

void somaVetores(double* v1, double escalar, double* v2, double* r, int n){ //v1 + escalar*v2
	int i;
	for(i = 0; i < n; i++){
		r[i] = v1[i] + escalar*v2[i];
	}
}

double* GC(double** A, double* x, double* b, int n){
	int i;
	for(i = 0; i < n; i++)
		x[i] = 0;
	
	double* r = malloc(sizeof(double)*n);
	
	double* Ax = malloc(sizeof(double)*n);
	double* Ar = Ax;
	double s;
	int ite = 0;
	do{
		multiplicaMatrizPorVetor(A,x,Ax,n); //Calcula Ax
		somaVetores(b,-1,Ax,r,n); //r = b + (-1*Ax)
		
		multiplicaMatrizPorVetor(A,r,Ar,n); //Calcula Ar
		s = multiplicaVetores(r,r,n) / multiplicaVetores(r,Ar,n); // s = r(t)*r / r(t)*Ar
		
		somaVetores(x,s,r,x,n); //x(k+1) = x(k) + sr
		
		
		ite++;
	}while(multiplicaVetores(r,r,n) > erroMax  &&  ite <= maxIte);
	printf("Terminou em %d iteracoes\n",ite);
	return x;
}


int main(){
	srand( 20162 );
	double* b = malloc(sizeof(double)*N);
	double** A = geraMatrizEscalonada(b);
	
	
	double* x = malloc(sizeof(double)*N);
	
	
	puts("----------Gradiente Conjugado-----------");
	double* y = GC(A,x,b,N);
	puts("----------------------------------------");
	puts("Gradiente Conjugado:");
	puts("x:");
	imprimeVetor(y);
	return 0;
}

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <likwid.h>
#include "immintrin.h"
#include "metodoGaus.c"


void calculaATAeATY(double* x, double* y,double** ATA, double* ATY, int p, int n){
	LIKWID_MARKER_START("antigo");
	int i,j;
	
	//Esse eh o numero máximo de elementos diferentes que eu tenho (ou o número de diagonais que existem na minha matriz)
	/*
	 Por exemplo:
	 
	 |S0 S1 S2|  Temos exatamente 4 diagonais, e 4 elementos diferentes
	 |S1 S2 S3|  Cada S pode ser calculado com: SN = Somatorio(x^N)
	 |S2 S3 S4|  Perceba que o número máximo de K para SK vai até 2*N - 1, no exemplo 2*3-1 = 5
	 
	*/ 
	
	
	
	ATA[0][0] = (double)p; //Este elemento tem sempre o valor de p
	
	//Estou usando este vetor para que em cada iteração ele contenha os valores de x^k  onde k é o número da iteração, inicialmente k = 1
	double* aux_vec_x = malloc(sizeof(double)*p);
	for(i = 0; i < p; i++){
		aux_vec_x[i] = x[i];
	}
	
	
	
	
	
	//Como o primeiro elmento de ATY depende do somatório de Xi ^ 0, já inicializamos
	ATY[0] = 0.0;
	for(i = 0; i < p; i++){
		ATY[0] += y[i];
	}
	
	
	
	double som;
	
	//Percorra as N primeira diagonais (até a metade) a partir de 1: No exemplo seria S1 e S2
	for(i = 1; i < n; ++i){ 
		som = 0.0;
		
		
		
		ATY[i] = 0.0;
		for(j = 0; j < p; ++j){
			som+=aux_vec_x[j]; //Calcula o somatorio de X^i pois seria futuramente o valor de Si
			
			ATY[i] += aux_vec_x[j]*y[j]; //Calcula o valor de ATY
			aux_vec_x[j]*=x[j]; //Calcula os valores de x^i para a próxima iteração
			
		}
		
		//Coloca o valor calculado numa diagonal
		/*
			Exemplo:
			Para i = 1; temos:
			ATA[1][0]
			ATA[0][1]
			 
			Para i = 2:
			ATA[2][0]
			ATA[1][1]
			ATA[0][2]
		*/
		for(j = 0; j <= i; ++j){
			ATA[i-j][j] = som;
		}
	}
	
	
	//Apartir daqui, já preenchemos a matriz até a metade (no exemplo até S2) e também os valores de ATY
	
	//Precisamos preencher o resto da matriz, ou seja, a partir da diagonal N+1 até N*2-1
	
	
	int maxN;
	
	for(i = 1; i < n; ++i){
		som = 0.0;
		
		//Faz a mesma coisa do laço anterior, porém não precisamos mais calcular ATY
		for(j = 0; j < p; ++j){
			som+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
		}
		
		maxN = n-i; //Na diagonal desta iteração temos maxN elementos
		//Similar ao laço das N primeiras diagonais tempos aqui uma maneira de percorrer similar:
		/*
			Para i = 0:
			ATA[1][2]
			ATA[2][1]
			*No exemplo essas são as posições de S3 
			
		*/
		
		for(j = 0; j < maxN; ++j){
			ATA[i+j][n-1-j] = som;
		}
		
	}
	
	LIKWID_MARKER_STOP("antigo");
}



void calculaATAeATY_otimizado(double* x, double* y,double** ATA, double* ATY, int p, int n){
	LIKWID_MARKER_START("novo");
	int i,j;
	

	
	
	ATA[0][0] = (double)p; 
	
	double* aux_vec_x = malloc(sizeof(double)*p);
	ATY[0] = 0.0;
	for(i = 0; i < p; i++){
		aux_vec_x[i] = x[i];
		ATY[0] += y[i];
	}
	
	
	
	double som, som2, som3, som4, som5, som6, som7, som8;
	
	for(i = 1; i+7 < n; i+=8){ 
		som  = 0.0;
		som2 = 0.0;
		som3 = 0.0;
		som4 = 0.0;
		som5 = 0.0;
		som6 = 0.0;
		som7 = 0.0;
		som8 = 0.0;
		
		ATY[i] = 0.0;
		ATY[i+1] = 0.0;
		ATY[i+2] = 0.0;
		ATY[i+3] = 0.0;
		ATY[i+4] = 0.0;
		ATY[i+5] = 0.0;
		ATY[i+6] = 0.0;
		ATY[i+7] = 0.0;
		for(j = 0; j < p; ++j){
			som+=aux_vec_x[j];
			ATY[i] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
			
			som2+=aux_vec_x[j];
			ATY[i+1] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
			
			som3+=aux_vec_x[j];
			ATY[i+2] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
			
			som4+=aux_vec_x[j];
			ATY[i+3] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
			
			som5+=aux_vec_x[j];
			ATY[i+4] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
			
			som6+=aux_vec_x[j];
			ATY[i+5] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
			
			som7+=aux_vec_x[j];
			ATY[i+6] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
			
			som8+=aux_vec_x[j];
			ATY[i+7] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
		}
		
		
		
		for(j = 0; j <= i; ++j){
			ATA[i-j][j] = som;
			ATA[i+1-j][j] = som2;
			ATA[i+2-j][j] = som3;
			ATA[i+3-j][j] = som4;
			ATA[i+4-j][j] = som5;
			ATA[i+5-j][j] = som6;
			ATA[i+6-j][j] = som7;
			ATA[i+7-j][j] = som8;
		}
		ATA[i+1-j][j] = som2;
		
		ATA[i+2-j][j] = som3;
		ATA[i+1-j][j+1] = som3;
		
		ATA[i+3-j][j] = som4;
		ATA[i+2-j][j+1] = som4;
		ATA[i+1-j][j+2] = som4;
		
		ATA[i+4-j][j] = som5;
		ATA[i+3-j][j+1] = som5;
		ATA[i+2-j][j+2] = som5;
		ATA[i+1-j][j+3] = som5;
		
		ATA[i+5-j][j] = som6;
		ATA[i+4-j][j+1] = som6;
		ATA[i+3-j][j+2] = som6;
		ATA[i+2-j][j+3] = som6;
		ATA[i+1-j][j+4] = som6;
		
		ATA[i+6-j][j] = som7;
		ATA[i+5-j][j+1] = som7;
		ATA[i+4-j][j+2] = som7;
		ATA[i+3-j][j+3] = som7;
		ATA[i+2-j][j+4] = som7;
		ATA[i+1-j][j+5] = som7;
		
		ATA[i+7-j][j] = som8;
		ATA[i+6-j][j+1] = som8;
		ATA[i+5-j][j+2] = som8;
		ATA[i+4-j][j+3] = som8;
		ATA[i+3-j][j+4] = som8;
		ATA[i+2-j][j+5] = som8;
		ATA[i+1-j][j+6] = som8;
	}
	for(i = i; i < n; ++i){ 
		som = 0.0;
		
		ATY[i] = 0.0;
		for(j = 0; j < p; ++j){
			som+=aux_vec_x[j]; 
			
			ATY[i] += aux_vec_x[j]*y[j];
			aux_vec_x[j]*=x[j];
		}
		
		
		
		for(j = 0; j <= i; ++j){
			ATA[i-j][j] = som;
		}
	}
	
	
	
	
	
	
	
	int maxN;
	
	for(i = 1; i+7 < n; i+=8){
		som  = 0.0;
		som2 = 0.0;
		som3 = 0.0;
		som4 = 0.0;
		som5 = 0.0;
		som6 = 0.0;
		som7 = 0.0;
		som8 = 0.0;
		
		for(j = 0; j < p; ++j){
			som+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
			
			som2+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
			
			som3+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
			
			som4+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
			
			som5+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
			
			som6+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
			
			som7+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
			
			som8+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
		}
		
		maxN = n-i-7;
		for(j = 0; j < maxN; ++j){
			ATA[i+j][n-1-j] = som;
			ATA[i+1+j][n-1-j] = som2;
			ATA[i+2+j][n-1-j] = som3;
			ATA[i+3+j][n-1-j] = som4;
			ATA[i+4+j][n-1-j] = som5;
			ATA[i+5+j][n-1-j] = som6;
			ATA[i+6+j][n-1-j] = som7;
			ATA[i+7+j][n-1-j] = som8;
		}
		ATA[i+j][n-1-j] = som;
		ATA[i+j+1][n-2-j] = som;
		ATA[i+j+2][n-3-j] = som;
		ATA[i+j+3][n-4-j] = som;
		ATA[i+j+4][n-5-j] = som;
		ATA[i+j+5][n-6-j] = som;
		ATA[i+j+6][n-7-j] = som;
		
		ATA[i+1+j][n-1-j] = som2;
		ATA[i+1+j+1][n-2-j] = som2;
		ATA[i+1+j+2][n-3-j] = som2;
		ATA[i+1+j+3][n-4-j] = som2;
		ATA[i+1+j+4][n-5-j] = som2;
		ATA[i+1+j+5][n-6-j] = som2;
		
		ATA[i+2+j][n-1-j] = som3;
		ATA[i+2+j+1][n-2-j] = som3;
		ATA[i+2+j+2][n-3-j] = som3;
		ATA[i+2+j+3][n-4-j] = som3;
		ATA[i+2+j+4][n-5-j] = som3;
		
		ATA[i+3+j][n-1-j] = som4;
		ATA[i+3+j+1][n-2-j] = som4;
		ATA[i+3+j+2][n-3-j] = som4;
		ATA[i+3+j+3][n-4-j] = som4;
		
		ATA[i+4+j][n-1-j] = som5;
		ATA[i+4+j+1][n-2-j] = som5;
		ATA[i+4+j+2][n-3-j] = som5;
		
		ATA[i+5+j][n-1-j] = som6;
		ATA[i+5+j+1][n-2-j] = som6;
		
		ATA[i+6+j][n-1-j] = som7;
		
	}
	for(i = i; i < n; ++i){
		som = 0.0;
		for(j = 0; j < p; ++j){
			som+=aux_vec_x[j];
			aux_vec_x[j]*=x[j];
		}
		
		maxN = n-i;
		for(j = 0; j < maxN; ++j){
			ATA[i+j][n-1-j] = som;
		}
		
	}
	
	LIKWID_MARKER_STOP("novo");
}

int main(){
	LIKWID_MARKER_INIT;
	int n;
	int p;
	
	
	scanf("%d %d", &n, &p);
	
	p = 1000; n = 50;
	
	double* x = malloc(sizeof(double)*p);
	double* y = malloc(sizeof(double)*p);
	
	
	
	
	
	int i;
	
	for(i = 0; i < p; i++){
		//scanf("%lf %lf", &x[i], &y[i]);
		x[i] = i+1;
		y[i] = (i+1)*(i+1) - 2*(i+1)+ 3;
	}
	
	
	double** A = malloc(sizeof(double*) * n);
	for(i = 0; i < n; i++){
		A[i] = malloc(sizeof(double)*n);
	}
	
	
	
	
	
	double* z = malloc(sizeof(double)*n);
	
	
	calculaATAeATY_otimizado(x,y,A,z,p,n);
	
	
	puts("ATA:\n");
	for(i = 0; i < 3; i++){
		printf("|%lf|\n",A[i][0]);
		printf("|%lf|\n",A[n-1][i]);
	}
	
	puts("\nATY:\n");
	imprimeVetor(z,n);
	
	puts("-----------------------------------------------");
	
	calculaATAeATY(x,y,A,z,p,n);

	puts("ATA:\n");
	for(i = 0; i < 3; i++){
		printf("|%lf|\n",A[i][0]);
		printf("|%lf|\n",A[n-1][i]);
	}
	
	puts("\nATY:\n");
	imprimeVetor(z,n);
	
	
	
	
	
	
	double* reta = malloc(sizeof(double)*n);
	
	for(i = 0; i < n; i++){
		reta[i] = 0.0;
	}
	
	metodoDeGauss(A,z,n);
	
	retroSubstituicao(A,reta,z,n);
	
	puts("\n\nResultado:\n");
	imprimeVetor(reta,n);
	
	LIKWID_MARKER_CLOSE;
	
	return 0;
	
}

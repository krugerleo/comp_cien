#include <stdlib.h>
#include <stdio.h>


#define MAX 100

double* tabela_dif(double* vec_x, double* vec_y, int N){
	int i;
	int j;
	
	double* D = malloc(sizeof(double) * N);
	//fprintf(stderr,"%d\n",N);
	
	for(i = 0; i < N; ++i){
		D[i] = vec_y[i];
	}
	
	for(i = 1; i < N; ++i){
		for(j = 0; j < N-i; ++j){
			D[i*N+j] =   (D[(i-1)*N+j+1] - D[(i-1)*N+j])/(vec_x[j+i]-vec_x[j]);
		}
	}
	return D;
	
	
}

double p(double x, double* vec_x, double* D, int N){
	int i,j;
	double ret = D[0];
	double aux;
	for(i = 1; i < N; ++i){
		aux = D[i*N];
		//printf("p(n) = %lf",aux);
		for(j = 0; j < i; ++j){
			aux*=(x-vec_x[j]);
			//printf(" * %lf",aux);
		}
		//printf("\n");
		ret+= aux;
	}
	
	return ret;
	
}


void imprimeD(double* D, int N){
	int i, j;
	
	for(j = 0; j < N; ++j){
		for(i = 0; i < N-j; ++i){
			printf("%lf ",D[i*N+j]);
		}
		printf("\n");
	}
	
	
}


int main(){
	double* x = malloc(sizeof(double)*MAX);
	double* y = malloc(sizeof(double)*MAX);

	int N;
	
	scanf("%d",&N);
	
	int i;
	for(i = 0; i < N; i++){
		scanf("%lf %lf",&x[i],&y[i]);
	}
	
	/*
	for(i = 0; i < N; i++){
		printf("p%d : (%lf,%lf)\n",i,x[i],y[i]);
	}
	*/
	
	double* D = tabela_dif(x,y,N);
	
	imprimeD(D,N);
	
	double ask;
	scanf("%lf",&ask);
	
	printf("f(%lf) = %lf \n", ask,p(ask,x,D,N));
	
	
	
}

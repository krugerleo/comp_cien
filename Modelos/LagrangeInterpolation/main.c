#include <stdlib.h>
#include <stdio.h>


#define MAX 100

double Lag_i(int i, int N, double x, double* vec_x){
	int j;
	double ret = 1;
	
	for(j = 0; j < i; ++j){
		ret *= (x-vec_x[j])/(vec_x[i]-vec_x[j]);
	}
	
	for(j = i+1; j < N; ++j){
		ret *= (x-vec_x[j])/(vec_x[i]-vec_x[j]);
	}
	
	return ret;
}

double p(double x, double* vec_x, double* vec_y, int N){
	int i;
	double ret = 0;
	for(i = 0; i < N; i++){
		ret+= Lag_i(i,N,x,vec_x)*vec_y[i];
	}
	return ret;
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
	
	double ask;
	scanf("%lf",&ask);
	
	printf("f(%lf) = %lf \n", ask,p(ask,x,y,N));
}

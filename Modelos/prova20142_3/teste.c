#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <likwid.h>
#include "immintrin.h"

//#define N 200000

int main(int argc, char *argv[]){
	
	
	
	LIKWID_MARKER_INIT;
	
	int i,j;
	const int N=100000000;
	double val;
	srand(1); // semente random
	
	double* A = malloc(sizeof(double)*N);
	double* B = malloc(sizeof(double)*N);
	double* C = malloc(sizeof(double)*N);
	double* D = malloc(sizeof(double)*N);
	/* */
	
	LIKWID_MARKER_START("lento");
	for(i = 0; i < N; ++i){
		A[i] = B[i] + C[i] * D[i];
	} 
	
	
	LIKWID_MARKER_STOP("lento");
	
	 /* */
	
	
	
	/* */
	
	double* E = aligned_alloc(32,sizeof(double)*N);
	double* F = aligned_alloc(32,sizeof(double)*N);
	double* G = aligned_alloc(32,sizeof(double)*N);
	double* H = aligned_alloc(32,sizeof(double)*N);
	
	
	LIKWID_MARKER_START("rapido");
	__m256d aux1;
	__m256d aux2;
	__m256d aux3;
	__m256d result;
	
	for(i = 0; i+3 < N; i+=4){
		aux1 = _mm256_setr_pd(F[i],F[i+1],F[i+2],F[i+3]);
		aux1 = _mm256_setr_pd(G[i],G[i+1],G[i+2],G[i+3]);
		aux1 = _mm256_setr_pd(H[i],H[i+1],H[i+2],H[i+3]);
		
		result = _mm256_fmadd_pd(aux2,aux3,aux1);
		E[i]   = ((double*)&result)[0];
		E[i+1] = ((double*)&result)[1];
		E[i+2] = ((double*)&result)[2];
		E[i+3] = ((double*)&result)[3];
	}
	for(i = i; i < N; ++i){
		E[i] = F[i] + G[i] * H[i];
	}

	LIKWID_MARKER_STOP("rapido");
	
	 /* */
	
	
	
	
	LIKWID_MARKER_CLOSE;
	
	return 0;
}

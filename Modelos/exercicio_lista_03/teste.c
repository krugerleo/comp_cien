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
	const int N=350;
	double val;
	srand(1); // semente random
	
	
	/* */
	LIKWID_MARKER_START("lento");
	
	
	double mat[N][N],s[N][N];
	// ...
	
	for(i=0; i<N ; ++i) {
		val = (double)rand(); // numero random
		for(j=0; j<N; ++j) {
			mat[j][i] = (s[j][i] / val) - cos(j%256) * cos(j%256);
		}
	}
	LIKWID_MARKER_STOP("lento");
	
	fprintf(stderr,">>>%lf \n>>>%lf\n", mat[20][10],mat[N-1][N-1]);
	 /* */
	
	
	
	/* */
	LIKWID_MARKER_START("rapido");
	
	double* cos2 = malloc(sizeof(double)*256);
	double auxCos;
	
	for(i = 0; i < 256; i++){
		auxCos = cos(i);
		cos2[i] = auxCos*auxCos;
	}
	
	
	double* mat2 = malloc(sizeof(double)*N*N);
	double* s2 = malloc(sizeof(double)*N*N);
	
	for(i = 0; i+15 < N; i+=16){
		val = (double)rand(); // numero random
		for(j = 0; j < N; ++j){
			mat2[j*N+i] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+1] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+2] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+3] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+4] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+5] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+6] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+7] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+8] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+9] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+10] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+11] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+12] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+13] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+14] = (s2[j*N+i] / val) - cos2[j%256];
			mat2[j*N+i+15] = (s2[j*N+i] / val) - cos2[j%256];
			
		}
	}
	for(i = i; i < N; ++i){
		val = (double)rand(); // numero random
		for(j = 0; j < N; ++j){
			mat2[j*N+i] = (s2[j*N+i] / val) - cos2[j%256];
		}
	}

	LIKWID_MARKER_STOP("rapido");
	
	free(cos2);
	free(s2);
	free(mat2);
	
	fprintf(stderr,">>>%lf \n>>>%lf\n",mat2[20*N+10],mat2[(N-1)*N + N - 1]);
	 /* */
	
	
	
	
	LIKWID_MARKER_CLOSE;
	
	return 0;
}

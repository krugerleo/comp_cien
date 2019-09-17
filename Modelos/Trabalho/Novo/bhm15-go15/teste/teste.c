#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <likwid.h>

#define N 200000

int main(){
	
	
	
	LIKWID_MARKER_INIT;
	
	int i;
	double a[N], b[N],c[N], d[N];
	
	LIKWID_MARKER_START("Marcador 1");
	
	
	for(i = 0; i < N; ++i){
		a[i] = b[i] + c[i] * d[i];
	}

	LIKWID_MARKER_STOP("Marcador 1");
	
	LIKWID_MARKER_START("Marcador 2");
	
	for(i = 0; i < N; i+=2){
		a[i] = b[i] + c[i] * d[i];
		//a[i+1] = b[i+1] + c[i+1] * d[i+1];
	}

	LIKWID_MARKER_STOP("Marcador 2");
	
	LIKWID_MARKER_CLOSE;
}

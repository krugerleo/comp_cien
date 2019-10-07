#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "method.h"

int gaussSeidel (sys_t *S, double time[], unsigned int IT ) {
	int i, j, k;
	double hy = PI / (S->ny + 1);
	double hx = PI / (S->nx + 1);

	for (k = 1; k < IT; k++) {
		time[k] = timestamp();
		for (j = 1; j <= S->ny; j++) {
			for(i = 1; i <= S->nx; i++)
			S->U[j*S->nx+2+i] = (S->U[j*S->nx+2+i+1]*coefum(hx,hy) 
					+ S->U[j*S->nx+2+i-1]*coefdois(hx,hy)
					+ S->U[j+1*S->nx+2+i]*coeftres(hx,hy)
					+ S->U[j-1*S->nx+2+i+1]*coefqua(hx,hy)
					+ calFunc(i,j)) / divi(hx,hy);
		}
		time[k] = timestamp() - time[k];
	}
	return k;
}

double calFunc(int x, int y){
	return (double)(4*pow(PI,2)*(sin(2*PI*(double)x)*sinh(PI*(double)y)+sin(2*PI*(PI-(double)x))*sinh(PI*(PI-(double)y))));
}

sys_t* inic(unsigned int nx, unsigned int ny){
	sys_t *S = (sys_t *) malloc(sizeof(sys_t));
	S->nx = nx;
	S->ny = ny;

	for(unsigned int i = 0; i < ny+2; i++){
		S->U[i*nx+2] = 0;
		S->U[i*nx+2+nx+1] = 0;
	}

	for(unsigned int i = 0; i < nx+2; i++){
		S->U[i] = sin(2.0*PI*(double)i)*sinh(pow(PI,2));
		S->U[ny+2+i] = sin(2.0*PI*(PI-(double)i))*sinh(pow(PI,2));
	}
	return (S);
}

double coefum(double hx, double hy){
	return (double)(hx*pow(hy,2) - (2.0*pow(hy,2)));
}

double coefdois(double hx, double hy){
	return (double)(-1.0*hx-pow(hy,2)-(2.0*pow(hy,2)));
}

double coeftres(double hx, double hy){
	return (double)(pow(hx,2)*hy - (2.0*pow(hx,2)));
}

double coefqua(double hx, double hy){
	return (double)(-1.0*pow(hx,2)*hy-(2.0*pow(hy,2)));
}

double divi(double hx, double hy){
	return (double)(2.0*pow(hx,2)-2.0*pow(hy,2)-2.0*pow(hx,2)*pow(hy,2)*ENE);
}	

void prnSistLinear (sys_t *SL)
{
  for(int i=0; i < SL->ny; ++i) {
    printf("\n");
    for(int j=0; j < SL->nx; ++j)
      printf ("%10g", SL->U[i*SL->nx+2+j]);
  }
  printf("\n\n");
}

double timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

void inicdois(double *U,int nx,int ny){
	for(i=0; i<nx+2;i++){
		U[i][0] = 1.1;
		U[i][ny+1] = 1.1;
	}

	for
}
int main (int argc, char *argv[]) {
    if (argc != 9) {
        printf("Use: -nx <Num pontos em X> -ny <Num potnso em y> -i <Num max de iterações> -o <Caminho para o arquivo de saida>\n");
        exit(-1);
    }
    for(int i=0; i<argc;i++){
    	printf("%s\n",argv[i] );
    }


    int NX = atoi(argv[2]);
    int NY = atoi(argv[4]);
    int IT = atol(argv[6]);
    sys_t *S;
    double U[NX+2][NY+2];
    inicdois(U,NX,NY);
    printf("%d e %d\n",NX,NY );
    S = inic(NX,NY);
    prnSistLinear(S);

    exit(0);
}

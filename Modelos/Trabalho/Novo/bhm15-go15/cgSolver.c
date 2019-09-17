#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <likwid.h>

#include "immintrin.h"

//#define checaDivisaoPorZero

double *antA = NULL;

int maxIte;
double tolerancia = -1;

FILE *output;

double FLT_EPSILON;



typedef struct _tempo{
	double min;
	double max;
	double somatorio;
} tempo;

void atualizaEstruturaTempo(double novoValor, tempo* estrutura_tempo){

	//Atualiza somatorio (usado para calcular tempo medio)
		estrutura_tempo->somatorio += novoValor;
	//

	//Atualiza tempo minimo
		if(estrutura_tempo->min == -1.0)
			estrutura_tempo->min = novoValor;
		else if(novoValor < estrutura_tempo->min){
			estrutura_tempo->min = novoValor;
		}
	//

	//Atualiza tempo maximo
		if(estrutura_tempo->max == -1.0)
			estrutura_tempo->max = novoValor;
		else if(novoValor > estrutura_tempo->max){
			estrutura_tempo->max = novoValor;
		}
	//
}

//double staticX[N] = {10000.2 ,20000.3 ,30000.4};


/***********************
 * N: tamanho do sistema linear
 * k: numero da diagonal, 0 = diagonal principal, 1 = acima/abaixo da diagonal, 2 = ...
 * nBandas: numero de bandas do sistema linear
 * diag: vetor para armazenar os valores da diagonal. Deve ser alocado por quem chama a função.
 ***********************/
int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int nBandas, double *diag )
{
  int i;

  if ( !diag || N < 3 || nBandas > N/2 || k < 0 || k > nBandas )
    return (-1);

  /* garante valor dominante para diagonal principal */
  double fator = (k == 0) ? ((double)(nBandas-1)) : (0.0);

  double invRandMax = 1.0 / (double)RAND_MAX;

  for (i=0; i < N-k; ++i)
  {
    diag[i] = fator + (double)rand() * invRandMax;
  }

  return (0);
}

#define SNN(i) i*(i+1)/2

#define blockPadding 32

#define A(i,j) A[(i*(n+32)  + j)]


double* geraMatriz(double* b, int bandas, int n){

	if(!(bandas % 2)){
		fprintf(stderr,"Número de bandas não é ímpar\n");
		exit(-1);
	}
	
	double* A = NULL;
	
	int linhasMatriz = (bandas/2) + 1;
	
	
	//Se preocupa com o alinhamento (endereco do vetor deve ficar num multiplo de 64)
	A = aligned_alloc(32,sizeof(double)*linhasMatriz*(n+blockPadding));
	
	//
	
	
	int i;
	for(i = 0; i < linhasMatriz; ++i)
	{
		if(generateRandomDiagonal(n,i, bandas, ((double*)&A[i*(n+blockPadding)]))){
			fprintf(stderr,"Erro ao gerar matriz\n");
			exit(-1);
		}
		
	}

	return A;
}

//f(x) = 4π² ( sin(2πx) + sin(2π(π-x)) )
double funcaoDeB(double xi){
	double auxFloat1 = 2*M_PI*xi;
	double auxFloat2 = 2*M_PI*(M_PI-xi);

	return 4*M_PI*M_PI * ( sin(auxFloat1) + sin(auxFloat2) );
}




double* geraB(int n){
	int i;
	double* b = malloc(sizeof(double)*n);
	for(i = 0; i < n; i++)
		b[i] = funcaoDeB( (double) (i * M_PI / n) );// para todo i = 0:n-1,

	return b;
}

double timestamp(void){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}


void multiplicaMatrizDeBandaPorVetor(double* A, double* v, double* r, int bandas, int n, tempo* t){
	int i,j, ij;
	double *p;
	double aux_tempo = timestamp();
	int linhasMatriz = (bandas>>1) + 1;
	
	int k;
	
	__m256d auxVec1;
	__m256d auxVec2;
	__m256d auxVec3;
	__m256d result;
	
	LIKWID_MARKER_START("Multiplicação Matriz de Banda X Vetor");
	p = A;
	
	for(i = 0; i+4 < n; i+=4){
		auxVec1 = _mm256_set_pd(p[i],p[i+1],p[i+2],p[i+3]);
		auxVec2  = _mm256_set_pd(v[i],v[i+1],v[i+2],v[i+3]);
		result = _mm256_mul_pd(auxVec1, auxVec2);
		r[i] = ((double*)&result)[3];
		r[i+1] = ((double*)&result)[2];
		r[i+2] = ((double*)&result)[1];
		r[i+3] = ((double*)&result)[0];
	}
	
	for(i = i; i < n; i++){
		r[i] = p[i]*v[i];
	}
	
	
	int l;
	int bitSee = 8176;
	for(l = 1; l+1 < linhasMatriz; l+=2){
		for(j = 0; j+linhasMatriz+6 < n; j+=8){
			
			for(i = l; i < l+2; ++i){
				p = &A(i,j);
				//printf("|	>>%d\n",(int)p & bitSee);
				auxVec1 = _mm256_setr_pd(p[0],p[1],p[2],p[3]);
				
				ij= i+j;
				auxVec2  = _mm256_setr_pd(v[ij],v[ij+1],v[ij+2],v[ij+3]);
				auxVec3  = _mm256_setr_pd(r[j],r[j+1],r[j+2],r[j+3]);
				result = _mm256_fmadd_pd(auxVec1, auxVec2,auxVec3);
				r[j] = ((double*)&result)[0];
				r[j+1] = ((double*)&result)[1];
				r[j+2] = ((double*)&result)[2];
				r[j+3] = ((double*)&result)[3];
				
				//printf("|		>>%d\n",(int)&v[j] & bitSee);
				//printf("|			>>%d\n",(int)&r[ij] & bitSee);
				auxVec2  = _mm256_setr_pd(v[j],v[j+1],v[j+2],v[j+3]);
				auxVec3  = _mm256_setr_pd(r[ij],r[ij+1],r[ij+2],r[ij+3]);
				result = _mm256_fmadd_pd(auxVec1, auxVec2,auxVec3);
				r[ij] = ((double*)&result)[0];
				r[ij+1] = ((double*)&result)[1];
				r[ij+2] = ((double*)&result)[2];
				r[ij+3] = ((double*)&result)[3];
				
				j+=4;
				ij+=4;
				p = &A(i,j);
				auxVec1 = _mm256_setr_pd(p[0],p[1],p[2],p[3]);
				
				auxVec2  = _mm256_setr_pd(v[ij],v[ij+1],v[ij+2],v[ij+3]);
				auxVec3  = _mm256_setr_pd(r[j],r[j+1],r[j+2],r[j+3]);
				result = _mm256_fmadd_pd(auxVec1, auxVec2,auxVec3);
				r[j] = ((double*)&result)[0];
				r[j+1] = ((double*)&result)[1];
				r[j+2] = ((double*)&result)[2];
				r[j+3] = ((double*)&result)[3];
				
				
				auxVec2  = _mm256_setr_pd(v[j],v[j+1],v[j+2],v[j+3]);
				auxVec3  = _mm256_setr_pd(r[ij],r[ij+1],r[ij+2],r[ij+3]);
				result = _mm256_fmadd_pd(auxVec1, auxVec2,auxVec3);
				r[ij] = ((double*)&result)[0];
				r[ij+1] = ((double*)&result)[1];
				r[ij+2] = ((double*)&result)[2];
				r[ij+3] = ((double*)&result)[3];
				j-=4;
			}
		
		}
	}
	
	
	for(j = 0; j+linhasMatriz+6 < n; j+=8){
		
		for(i = l; i < linhasMatriz; ++i){
			p = &A(i,j);
			auxVec1 = _mm256_setr_pd(p[0],p[1],p[2],p[3]);
			
			ij= i+j;
			auxVec2  = _mm256_setr_pd(v[ij],v[ij+1],v[ij+2],v[ij+3]);
			auxVec3  = _mm256_setr_pd(r[j],r[j+1],r[j+2],r[j+3]);
			result = _mm256_fmadd_pd(auxVec1, auxVec2,auxVec3);
			r[j] = ((double*)&result)[0];
			r[j+1] = ((double*)&result)[1];
			r[j+2] = ((double*)&result)[2];
			r[j+3] = ((double*)&result)[3];
			
			
			auxVec2  = _mm256_setr_pd(v[j],v[j+1],v[j+2],v[j+3]);
			auxVec3  = _mm256_setr_pd(r[ij],r[ij+1],r[ij+2],r[ij+3]);
			result = _mm256_fmadd_pd(auxVec1, auxVec2,auxVec3);
			r[ij] = ((double*)&result)[0];
			r[ij+1] = ((double*)&result)[1];
			r[ij+2] = ((double*)&result)[2];
			r[ij+3] = ((double*)&result)[3];
			
			j+=4;
			ij+=4;
			p = &A(i,j);
			auxVec1 = _mm256_setr_pd(p[0],p[1],p[2],p[3]);
			
			auxVec2  = _mm256_setr_pd(v[ij],v[ij+1],v[ij+2],v[ij+3]);
			auxVec3  = _mm256_setr_pd(r[j],r[j+1],r[j+2],r[j+3]);
			result = _mm256_fmadd_pd(auxVec1, auxVec2,auxVec3);
			r[j] = ((double*)&result)[0];
			r[j+1] = ((double*)&result)[1];
			r[j+2] = ((double*)&result)[2];
			r[j+3] = ((double*)&result)[3];
			
			
			auxVec2  = _mm256_setr_pd(v[j],v[j+1],v[j+2],v[j+3]);
			auxVec3  = _mm256_setr_pd(r[ij],r[ij+1],r[ij+2],r[ij+3]);
			result = _mm256_fmadd_pd(auxVec1, auxVec2,auxVec3);
			r[ij] = ((double*)&result)[0];
			r[ij+1] = ((double*)&result)[1];
			r[ij+2] = ((double*)&result)[2];
			r[ij+3] = ((double*)&result)[3];
			j-=4;
		}
	
	}
		
	k = j;	
	for(i = 1; i < linhasMatriz; ++i){
		for(j = k; i+j < n; ++j){
			r[j]+=A(i,j)*v[i+j];
			r[i+j]+=A(i,j)*v[j];
		}
	}
	
	atualizaEstruturaTempo(timestamp() - aux_tempo, t);
	LIKWID_MARKER_STOP("Multiplicação Matriz de Banda X Vetor"); 
}


double multiplicaVetores(double* v1, double* v2, int n, tempo* t){
	
	double soma = 0.0;
	
	__m256d aux1;
	__m256d aux2;
	__m256d resultado;
	
	int i;
	double aux_tempo = timestamp();
        
	LIKWID_MARKER_START("Multiplicação Vetor X Vetor");
	
	
	for(i = 0; i+3 < n; i += 4){
		aux1 = _mm256_set_pd(v1[i],v1[i+1],v1[i+2],v1[i+3]);
		aux2  = _mm256_set_pd(v2[i],v2[i+1],v2[i+2],v2[i+3]);
		
		
		resultado = _mm256_mul_pd(aux1, aux2);
		soma  += ((double*)&resultado)[0] + ((double*)&resultado)[1] + ((double*)&resultado)[2] + ((double*)&resultado)[3];
	}
	
	
	for(i = i; i < n; i++){
		soma+=v1[i]*v2[i];
	}
	
	atualizaEstruturaTempo(timestamp() - aux_tempo, t);

    LIKWID_MARKER_STOP("Multiplicação Vetor X Vetor");

	return soma;
}

void somaVetores(double* v1, double escalar, double* v2, double* r, int n){ //r = v1 + escalar*v2
	int i;
	for(i = 0; i < n; i++){
		r[i] = v1[i] + escalar*v2[i];
	}
}

void GC(double* A, double* x, double* b, int bandas, int n){
	int i;
	double auxTempo_gc;
	double auxTempo_r;


	tempo tempo_gc;
	tempo_gc.min = -1.0;
	tempo_gc.max = -1.0;
	tempo_gc.somatorio = 0.0;

	tempo tempo_r;
	tempo_r.min = -1.0;
	tempo_r.max = -1.0;
	tempo_r.somatorio = 0.0;

	tempo tempo_mmv;
	tempo_mmv.min = -1.0;
	tempo_mmv.max = -1.0;
	tempo_mmv.somatorio = 0.0;

	tempo tempo_mvv;
	tempo_mvv.min = -1.0;
	tempo_mvv.max = -1.0;
	tempo_mvv.somatorio = 0.0;

	//x0 = 0
	for(i = 0; i < n; i++)
		x[i] = 0;

	//v = b
	double* v = malloc(sizeof(double)*n);

	//r = b
	double* r = malloc(sizeof(double)*n);

	//como, inicialmente, r = b = v
	for(i = 0; i < n; i++){
		v[i] = b[i];
		r[i] = b[i];
	}

	//aux = r^t * r
	double aux = multiplicaVetores(r,r,n,&tempo_mvv);

	double* z = malloc(sizeof(double)*n);
	double s;
	double aux1;
	double m;


	double* residuo = malloc(sizeof(double)*n);

	double* norma = malloc(sizeof(double)*maxIte+1);
	norma[0] = sqrt(aux); //Norma de b (erro inicial pois x = 0);

	double erroAproximado;







	void encerraFunc(){
			int j;

			fputs("###########\n",output);

			fprintf(output,"# Tempo Método CG: %.14g %.14g %.14g \n", tempo_gc.min, (tempo_gc.somatorio/(double)(i)), tempo_gc.max);
			fprintf(output,"# Tempo Resíduo: %.14g %.14g %.14g \n#\n", tempo_r.min, (tempo_r.somatorio/(double)(i)), tempo_r.max);

			fprintf(output,"# Tempo MMV: %.14g %.14g %.14g \n", tempo_mmv.min, (tempo_mmv.somatorio/(double)(i*2)), tempo_mmv.max);
			fprintf(output,"# Tempo MVV: %.14g %.14g %.14g \n#\n", tempo_mvv.min, (tempo_mvv.somatorio/(double)((i*4)+1)), tempo_mvv.max);

			fprintf(output,"# Norma Euclidiana do Residuo e Erro aproximado\n");
			for(j = 1; j <= i; ++j) fprintf(output,"#  i=%d: %.14g\n", j, norma[j]); // i = 0 contem a norma da primeira iteracao ( ||b|| )

			fputs("###########\n",output);

			//Libera os vetores
			free(v);
			free(r);
			free(z);
			free(residuo);
			free(norma);
	}

	i = 0;
	do{
     	/*Calculo de X*/
		++i;

		auxTempo_gc = timestamp();
		
		multiplicaMatrizDeBandaPorVetor(A,v,z,bandas,n,&tempo_mmv); //z = Av

		#ifdef checaDivisaoPorZero
		if(FLT_EPSILON > fabs(multiplicaVetores(v,z,n,&tempo_mvv))) {
			fprintf(stderr,"Erro: divisão por zero\n");
			exit(-1);
		}
		#endif
		s = aux / multiplicaVetores(v,z,n,&tempo_mvv); // s = aux / v^t * z

		somaVetores(x,s,v,x,n); // x = x + s*v

		somaVetores(r,-s,z,r,n); //r = r - s*z

		aux1 = multiplicaVetores(r,r,n,&tempo_mvv); //aux1 = r^t * r


		/*Calculo do resíduo*/

		auxTempo_r = timestamp();

		multiplicaMatrizDeBandaPorVetor(A,x,residuo,bandas,n,&tempo_mmv); //residuo = A*x
		somaVetores(b,-1,residuo,residuo,n); //residuo = b - residuo;
		norma[i] = sqrt( multiplicaVetores(residuo,residuo,n,&tempo_mvv)); //norma euclidiana do residuo
		erroAproximado = fabs(norma[i-1] - norma[i]); //erro aproximado

		atualizaEstruturaTempo(timestamp() - auxTempo_r,&tempo_r);




		if( erroAproximado  < tolerancia ){
			atualizaEstruturaTempo(timestamp() - auxTempo_gc,&tempo_gc);
			encerraFunc();

			return; // Saiu quando o erro aproximado foi menor que a tolerancia
		}

		#ifdef checaDivisaoPorZero
		if(FLT_EPSILON > fabs(aux)) {
			fprintf(stderr,"Erro: divisão por zero\n");
			exit(-1);
		}
		#endif
		m = aux1/aux;
		aux = aux1;
		somaVetores(r,m,v,v,n); //v = r + m*v

		atualizaEstruturaTempo(timestamp() - auxTempo_gc,&tempo_gc);

	}while(i < maxIte);

	encerraFunc();


	return;// Saiu com o máximo de iterações
}

int main(int argc, char **argv){
	LIKWID_MARKER_INIT;

	FLT_EPSILON = 1.0;

	while(1.0 - FLT_EPSILON != 1.0)
		FLT_EPSILON /= 2;

	FLT_EPSILON *= 2;

	int N;
	int bandas;

	if( (!argv[1]) || (!argv[2]) ){
		fprintf(stderr,"Falta de parâmetro(s) obrigatório(s)\n");
		exit(-1);
	}

	if (atoi(argv[1]))
		N = atoi(argv[1]);
	else {
		fprintf(stderr,"Parâmetro errado: Dimensão da matriz inválida\n");
		exit(-1);
	}


	if(atoi(argv[2]))
		bandas = atoi(argv[2]);
	else {
		fprintf(stderr,"Parâmetro errado: Número de bandas inválido\n");
		exit(-1);
	}

	maxIte = N;
	output = stdout;

	int i;
	for(i = 3; i < argc; i++){
		if( !strcmp(argv[i], "-i") ){
			if (atoi(argv[i+1]))
				maxIte = atoi(argv[i+1]);
			else{
				fprintf(stderr,"Parâmetro errado: Número de iterações inválido\n");
				exit(-1);
			}

			i++;
		}

		if( !strcmp(argv[i], "-t") ){
			if (atof(argv[i+1]))
				tolerancia = atof(argv[i+1]);
			else{
				fprintf(stderr,"Parâmetro errado: Número de tolerância inválido\n");
				exit(-1);
			}
			i++;
		}

		if( !strcmp(argv[i], "-o") ){
			output = fopen( argv[i+1], "w");
			i++;
		}
	}


	srand( 20162 );


	double* b = geraB(N);
	double* A = geraMatriz(b,bandas,N);



	double* x = malloc(sizeof(double)*N);

	GC(A,x,b,bandas,N);

	fprintf(output,"%d\n", N);
	for(i=0; i < N; i++) fprintf(output,"%.14g ", x[i]);
	fputs("\n",output);



	free(x);
	free(b);
	
	free(A);


	if( (output != stdout) && output)
		fclose(output);

	LIKWID_MARKER_CLOSE;

	exit(0);

}

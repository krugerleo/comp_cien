#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <likwid.h>
#define checaDivisaoPorZero


int maxIte;
double tolerancia = -1;

FILE *output;

double FLT_EPSILON;



typedef struct _tempo{
	double min;
	double max;
	double somatorio;
} tempo;

double timestamp(void){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

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


double** geraMatriz(double* b, int bandas, int n){

	if(!(bandas % 2)){
		fprintf(stderr,"Número de bandas não é ímpar\n");
		exit(-1);
	}
	double** A = malloc(sizeof(double*)*bandas);



	int k = 0;
	int i,j;
	for(i = (int)(bandas/2); i < bandas; ++i)
	{
		A[i] = malloc(sizeof(double)*n);
		//Melhorar para preencher soh os primeiros
		for(j = 0; j < n; j++)
				A[i][j] = 0.0;

		if(generateRandomDiagonal(n,k, bandas, A[i])){
			fprintf(stderr,"Erro ao gerar matriz\n");
			exit(-1);
		}


		if(k != 0)
		{
			A[ (bandas/2) - k ] = malloc(sizeof(double)*n);

			//Vamos "empurrar" os elementos para a direita
			for(j = k; j < n; j++)
			{
				A[ (bandas/2) - k ][j] = A[i][j-k];
			}

		}

		++k;
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

void multiplicaMatrizPorVetor(double** A, double* v, double* r, int n){
	int i,j;
	for(i = 0; i < n; i++){
		r[i] = 0.0;
		for(j = 0; j < n; j++){
			r[i] += A[i][j]*v[j];
		}
	}
}

void multiplicaMatrizDeBandaPorVetor(double** A, double* v, double* r, int bandas, int n, tempo* t){
	int i,j;
	int k;
	double aux_tempo;
	LIKWID_MARKER_START("Multiplicação Matriz de Banda X Vetor");
	
	aux_tempo = timestamp();
	for(j = 0; j < n; j++){
		r[j] = 0.0;
		k = j - bandas/2;

		for(i = 0; i < bandas ; i++){
			if(( (k+i) >= 0 ) && ( (k+i) < n )){
				r[j] += A[i][j] * v[k+i];
			}

		}
	}
	atualizaEstruturaTempo(timestamp() - aux_tempo, t);

	LIKWID_MARKER_STOP("Multiplicação Matriz de Banda X Vetor");	
}


double multiplicaVetores(double* v1, double* v2, int n, tempo* t){
	double soma = 0.0;
	int i;
	double aux_tempo;
	LIKWID_MARKER_START("Multiplicação Vetor X Vetor");
	
	aux_tempo = timestamp();
	for(i = 0; i < n; i++){
		soma += (v1[i]*v2[i]);
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


void GC(double** A, double* x, double* b, int bandas, int n){
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
	double** A = geraMatriz(b,bandas,N);



	double* x = malloc(sizeof(double)*N);

	GC(A,x,b,bandas,N);

	fprintf(output,"%d\n", N);
	for(i=0; i < N; i++) fprintf(output,"%.14g ", x[i]);
	fputs("\n",output);



	free(x);
	free(b);
	for(i = 0; i < bandas; i++)
		free(A[i]);
	free(A);


	if( (output != stdout) && output)
		fclose(output);

	LIKWID_MARKER_CLOSE;

	exit(0);

}

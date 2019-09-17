double funcao(double*A,double*b,double*x, unsigned int N,double maxErr, int maxIter ){
   double resSq, lambda, numerador, denominador; 
   unsigned int i, j, iter; 
   double *res = (double*) malloc (N*sizeof(double)); 
   double *vetorTemp = (double*) malloc (N*sizeof(double)); 

   memset(x, 0, N * sizeof(double));
    for(resSq = 0.0, i=0; i < N; ++i) { 
       res[i] = 0.0; 
       for(j=0; j < N; ++j) { 
          res[i] += A[i*N+j] * x[j]; 
       } 
       res[i] = b[i] ­ res[i]; 
       resSq += res[i] * res[i];
    }



	for(iter=0; iter < maxIter && sqrt(resSq) > maxErr; ++iter) {
		numerador = 0.0;
		denominador = 0;
		resSq = 0.0;
		lambda = numerador / denominador;
		for(i=0; i+3 < N; ++i) { 
	       vetorTemp[i]   = 0.0; 
	       vetorTemp[i+1] = 0.0; 
	       vetorTemp[i+2] = 0.0; 
	       vetorTemp[i+3] = 0.0; 
	
			for(j=0; j < N; ++j) { 
		          vetorTemp[i]   += A[i*N+j] * res[j]; 
		          vetorTemp[i+1] += A[(i+1)*N+j] * res[j]; 
		          vetorTemp[i+2] += A[(i+2)*N+j] * res[j]; 
		          vetorTemp[i+3] += A[(i+3)*N+j] * res[j]; 
		    }
		
			numerador += res[i] * res[i];
			denominador += vetorTemp[i] * res[i];
			x[i] += lambda * res[i];
			
			numerador += res[i+1] * res[i+1];
			denominador += vetorTemp[i+1] * res[i+1];
			x[i+1] += lambda * res[i+1];

			numerador += res[i+2] * res[i+2];
			denominador += vetorTemp[i+2] * res[i+2];
			x[i+2] += lambda * res[i+2];

			numerador += res[i+3] * res[i+3];
			denominador += vetorTemp[i+3] * res[i+3];
			x[i+3] += lambda * res[i+3];
		}
		for(i=i; i < N; ++i) { 
	       vetorTemp[i] = 0.0; 
	
			for(j=0; j < N; ++j) { 
		          vetorTemp[i] += A[i*N+j] * res[j]; 
		    }
		
			numerador += res[i] * res[i];
			denominador += vetorTemp[i] * res[i];
			x[i] += lambda * res[i];
		}
	
	    

		for(i=0; i+3 < N; i+=4) { 
			res[i] = 0.0; 
			res[i+1] = 0.0; 
			res[i+2] = 0.0; 
			res[i+3] = 0.0; 
			for(j=0; j < N; ++j){ 
				res[i] += A[i*N+j] * x[j]; 
				res[i+1] += A[(i+1)*N+j] * x[j]; 
				res[i+2] += A[(i+2)*N+j] * x[j]; 
				res[i+3] += A[(i+3)*N+j] * x[j]; 
			}
			res[i] = b[i] ­- res[i]; 
			resSq += res[i] * res[i];
			
			res[i+1] = b[i+1] ­- res[i+1]; 
			resSq += res[i+1] * res[i+1];
			
			res[i+2] = b[i+2] ­- res[i+2]; 
			resSq += res[i+2] * res[i+2];
			
			res[i+3] = b[i+3] ­- res[i+3]; 
			resSq += res[i+3] * res[i+3];
	    }
		for(i=i; i < N; ++i) { 
	       res[i] = 0.0; 
			for(j=0; j < N; ++j){ 
				res[i] += A[i*N+j] * x[j]; 
			} 
	       res[i] = b[i] ­- res[i]; 
	       resSq += res[i] * res[i];
	    }
	
	} // for (iter)

}

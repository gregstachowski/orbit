#include <math.h>

int gauss(double *a, double *x, int N) {
	int i, j, k;
	int max;
	int M = N+1; 
	double tmp; 

	/* printf("\n\n",M); */

	for(i=0;i<N;i++) {
	/*	 for( j=0;j<N+1;j++) {
			printf(" % 1.7f ", *(a+j+i*M));		
		} */
		*(x+i) = 0; 
		/* printf(" \t\t% 1.7f", *(x+i)); 
		printf("\n");  */
	}
	
	for(i=0;i<N;i++) {			/* eliminacja */ 
  		max=i; 
  		for(j=i+1;j<N;j++) {
	 		if(fabs(*(a+i+j*M))>fabs(*(a+i+max*M))) {
				max=j; 
			}
 		}

		for(k=i;k<N+1;k++) { 		/* zamiana wierszy wartoœciami */
			tmp=(*(a+k+i*M)); 
		 	*(a+k+i*M)=(*(a+k+max*M)); 
		 	*(a+k+max*M)=tmp; 
		} 
 
		if(*(a+i+i*M)==0) {
  			return(0);  		/* Uk³ad sprzeczny! */
 		}
 
 	 	for(j=i+1;j<N;j++) {
		 	for(k=N;k>=i;k--) {	/* mno¿enie wiersza j przez wspó³czynnik "zeruj¹cy": */
				*(a+k+j*M)-=((*(a+k+i*M)*(*(a+i+j*M)))/(*(a+i+i*M))); 
  			} 
		}
	}

/*	printf("\n");
	for(i=0;i<N;i++) {
		for( j=0;j<N+1;j++) {
			printf(" % 1.7f ", *(a+j+i*M));
		}
		printf("\n");
	} 
*/
	#ifdef TEST 
  	printf("MACIERZ TRÓJK¥TNA\n"); 
  	for(i=0;i<N;i++) { 
		for(j=0;j<=N;j++) {
			printf(" % 1.7f ",*(a+j+i*M)); 
			printf("\n"); 
		}
	}
 	#endif 

	/* reduckja wsteczna */

	for(j=N-1;j>=0;j--) { 
  		tmp=0; 
  		for(k=j+1;k<N;k++) {
			tmp+=*(a+k+j*M)*(*(x+k));
		}
  		*(x+j)=(*(a+N+j*M)-tmp)/(*(a+j+j*M)); 
  	} 

	return(1);  /* wszystko w porz¹dku! */
} 
 

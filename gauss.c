#include <math.h>

#define N 3

int gauss(double a[N][N+1], double x[N]) {
	int i, j, k;
	int max; 
	double tmp; 
	double y[N+1];

	for(i=0;i<N;i++) {
		for( j=0;j<N+1;j++) {
		/*	printf(" % 1.7f ", a[i][j]);*/
		y[j] = 0;
		}
		x[i] = 0;
		/* printf(" \t\t% 1.7f", x[i]); 
		printf("\n"); */
	}

	for(i=0;i<N;i++) {		/* eliminacja */ 
  		max=i; 
  		for( j=i+1;j<N;j++) {
	 		if(fabs(a[j][i])>fabs(a[max][i])) {
				max=j; 
			}
 		}

		for(k=i;k<N+1;k++) { 		/* zamiana wierszy wartoœciami */
			tmp=a[i][k]; 
		 	a[i][k]=a[max][k]; 
		 	a[max][k]=tmp; 
		} 
 
		if(a[i][i]==0) {
  			return(0);  		/* Uk³ad sprzeczny! */
 		}
 
 	 	for(j=i+1;j<N;j++) {
		 	for(k=N;k>=i;k--) {	/* mno¿enie wiersza j przez wspó³czynnik "zeruj¹cy": */
				a[j][k]-=((a[i][k]*a[j][i])/a[i][i]); 
  			} 
		}
 	}

/*	printf("\n");
	for(i=0;i<N;i++) {
		for( j=0;j<N+1;j++) {
			printf(" % 1.7f ", a[i][j]);
		}
		printf("\n");
	} */


	#ifdef TEST 
  	printf("MACIERZ TRÓJK¥TNA\n"); 
  	for(i=0;i<N;i++) { 
		for(j=0;j<=N;j++) {
			printf(" % 1.7f ",a[i][j]); 
			printf("\n"); 
		}
	}
 	#endif 

	/* reduckja wsteczna */

	for(j=N-1;j>=0;j--) { 
  		tmp=0; 
  		for(k=j+1;k<=N;k++) {
			/* printf("%d %d % 1.7f\n",j,k,y[k]); */
			tmp+=a[j][k]*y[k];
		}
  		y[j]=(a[j][N]-tmp)/a[j][j]; 
  	} 

	/* printf("\n"); */
	for(j=0;j<N;j++) {
		/* printf(" % 1.7f ",y[j]); */
		x[j] = y[j];
	}
	/* printf("\n"); */

	return(1);  /* wszystko w porz¹dku! */
} 
 

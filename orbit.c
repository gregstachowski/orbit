/* Orbit determination from  three  sets of observations (t, R.A., dec.) 
   assuming heliocentric orbit and geocentric observer (i.e. requires R)
   and using Laplacian methods. See "Astrodynamics: Orbit Determination,
   Space Navigation,  Celestial Mechanics,  Volume 1" by Samuel Herrick,
   Van Nostrand Reinhold, 1971;  esp.  Ch. 10 and 12 (page numbers refer 
   to the appropriate pages in this edition).

   Requires an ASCII datafile, eight (8) lines:

      comment
      header
      t(1)    R.A.(1)	dec(1)	Rx(1)	Ry(1)	Rz(1)
      t(1)    R.A.(1)   dec(1)  Rx(2)   Ry(2)   Rz(2)
      t(1)    R.A.(1)   dec(1)  Rx(3)   Ry(3)   Rz(3)
      k*
      alpha
      epsilon

   where t is in JD and R.A., dec. are both in decimal degrees,

   and k* = k = 1/(unit of canonical time) :

	k = 0.017 202 09895 if canonical time in days
        k = 0.001 239 444   if canonical time in seconds

   note: tau = k*(t - to).

   and alpha is the light-time constant:

	alpha = 0.0057755 days/au or alpha = 0.021275 sec/gu

   epsilon is the obliquity (inclination) of the ecliptic plane, 
   in decimal degrees.

   GSS, 12 Mar 1998, 3 Jun 1998, 7-9 Oct 1998 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* --- External functions --- */

extern double rad(double val);
extern int gauss(double*, double*, int N);
extern double estimrho(	double*, 
			double R_cos_psi, 
			double R_dot_R, 
			double A_prime, 
			double B_prime, 
			double E_);

/* --- Function definitions --- */

void trig(	double a[3], 
		double d[3], 
		double cos_a[3], 
		double cos_d[3],
		double sin_a[3],
		double sin_d[3]	) {

	int i;
	for (i=0;i<3;i++) {
		cos_a[i] = cos(rad(a[i]));
                cos_d[i] = cos(rad(d[i]));
                sin_a[i] = sin(rad(a[i]));
                sin_d[i] = sin(rad(d[i]));
	}
}

double determinant( 	double a[3],
			double b[3],
			double c[3] ) {

	/* works out:

	             | a1 b1 c1 |
	             | a2 b2 c2	|
	             | a3 b3 c3 |
	*/

	double det;

	det = (a[0]*b[1]*c[2]) + (a[1]*b[2]*c[0]) + (a[2]*b[0]*c[1]);
	det = det - (c[0]*b[1]*a[2]) - (c[1]*b[2]*a[0]) - (c[2]*b[0]*a[1]);	

	return(det);
}

double dot(double x[3], double y[3]) {

	/* Scalar product of two vectors, x.y  */
	/* Iloczyn skalarny dwoch wektorow, x.y */

	int i;
	double d;

	d = 0;
	for(i=0;i<3;i++) {
		d+=x[i]*y[i];
	}
	return(d);
}

void lad(       double cos_a[3], 
                double cos_d[3],
                double sin_a[3],
                double sin_d[3],
		double L[3][3],
		double A[3][3],
		double D[3][3]		) {

	int i;

	/* p.318 */

	for (i=0;i<3;i++) {
		L[i][0] = cos_d[i]*cos_a[i];
		L[i][1] = cos_d[i]*sin_a[i];
		L[i][2] = sin_d[i];

		A[i][0] = -1.0*sin_a[i];
		A[i][1] =      cos_a[i];
		A[i][2] =  0;

		D[i][0] = -1.0*sin_d[i]*cos_a[i];
		D[i][1] = -1.0*sin_d[i]*sin_a[i];
		D[i][2] =      cos_d[i];
	}

}

double kepler(double ecosE, double esinE, double ntau) {
	double fo, fp, E;
	
	printf("Solving Kepler ");
	E = ntau;
	printf("% 1.7f ",E);
	do {
		fo =  100*E + esinE*(100-100*cos(E)) - 100*ecosE*sin(E) - 100*ntau;
		fp =  100 + 100*esinE*sin(E) - 100*ecosE*cos(E);
		E += -1*fo/fp;
		E = fmod(E,2*M_PI);
		printf(">");
	} while(fabs(fo)>1e-6);
	printf("\n");
	return(E);
}

void readdata(	char filename[64], 
		char comment[128],
		char header[128],
		double t[3], 
		double a[3],
		double d[3],
		double R[3][3],
		double *k, 
		double *alpha,
		double *epsilon	   ) {

	FILE *input = NULL;
	int i;
	char s[128];

	input = fopen(filename, "r");
	fgets(comment, 128, input);	/* read one-line comment */
	fgets(header,  128, input);	/* read one-line header  */
	for (i=0;i<3;i++) {
		fgets(s, 128, input);
		sscanf(s, "%lf %lf %lf %lf %lf %lf", &t[i], &a[i], &d[i], \
		                           &R[i][0], &R[i][1], &R[i][2]);
	}
	fgets(s, 128, input);
	sscanf(s, "%lf", k);
	fgets(s, 128, input);
	sscanf(s, "%lf", alpha);
	fgets(s, 128, input);
	sscanf(s, "%lf", epsilon);
	fclose(input);
}

/* -------------------------------------- */

int main(int argc, char *argv[]) {

	char 	filename[64];
	char    comment[128];
	char	header[128];
        double 	t[3], a[3], d[3];
	double  tau3, tau1, tau13, denom;
	double  J[3], G[3];
	double 	cos_a[3], cos_d[3], sin_a[3], sin_d[3];
	double  L[3][3], A[3][3], D[3][3], R[3][3];
	double 	R_cos_psi, R_dot_R;
	double	L_dot[3], L_ddot[3], R_dot[3], R_ddot[3], L_dot_squared;
	double  k, alpha, epsilon;
	double	rho[3], rho_dot[3], rho_[3][3], r[3][3], r_dot[3][3], r_sq, r_cu; 
	double  cosd_dela[3], deld[3];

	double  ro = 0, ro_dot = 0, so_dot_sq;	/* see p. 169 [7B4-6] */

	double  u;	/* approximate solution for rho	   */
	double  E_, A_, B_, A_prime, B_prime;	/*	-- '' --		   */

	double  f[2] = {1, 1}, g[2] = {0, 0};	/* sums for f and g series 	   */
	double  f_[2][3];
	double  fcoeff[4], gcoeff[4];		/* coefficients for f and g series */

	double	gin1[3][4],  delr[3];		/* for Gaussian elimination 	   */
	double  gin2[4][5], delrd[4];
	double  delrdo[3];

	double  aj[2], bj[2]; /* p. 396 */

	double  coseps, sineps, dummy;
	double  V_[3], U_[3], Do, vo, uo, Uzeo, Vzeo, lo;
	double  E[2], e, M, p, incl;
	double u_;

	char    rstr[4] = "xyz";
	int	i, j, ncorr=0;

	/* Read datafile */

	if (argc != 2) {
		printf("\nsyntax: orbit <datafilename>\n\n");
		exit(1);
	}	
	else {
		sscanf(argv[1],"%s",filename);
	}
	
	printf("\nreading from: %s\n",filename);
	readdata(filename, comment, header, t, a, d, R, &k, &alpha, &epsilon);
	printf("\n%s\n", comment);
	printf("+--Datafile---------------------------------------+\n");
	printf("    %s", header);
	for (i=0;i<3;i++) {
		printf(" %d: %2.4f\t%3.4f\t\t% 03.4f\t% .4f\t% .4f\t% .4f\n", \
		             i, t[i], a[i], d[i], R[i][0], R[i][1], R[i][2]);
	}
	printf(" 4: %.11f\n", k);
	printf(" 5: %.11f\n", alpha);
	printf(" 6: %2.11f\n", epsilon);

	/* Calculate canonical time interval [8B17] */
	/* p.232 */

	tau1  = -1*k*(t[1] - t[0]); /* tau12 = -tau21 = -tau1 = k(t2 - t1) */
	tau3  =    k*(t[2] - t[1]); /* tau23 =           tau3 = k(t3 - t2) */
	tau13 =  tau3 - tau1;       /* tau13 = tau3 - tau1 */

	/* Calculate cos a, cos d, sin a, sin d */

	trig(a, d, cos_a, cos_d, sin_a, sin_d);
	printf("+--cos_a, cos_d, sin_a, sin_d---------------------+\n");
	for (i=0;i<3;i++) {
                printf(" %d: % 1.7f, % 1.7f", i, cos_a[i], cos_d[i]);
	        printf(" % 1.7f, % 1.7f\n",  sin_a[i], sin_d[i]);
        }
	
	/* Calculate L, A, D */

	lad(cos_a, cos_d, sin_a, sin_d, L, A, D);

        printf("+--L1 2 3-----------------------------------------+\n");
        for (j=0;j<3;j++) {
		printf(" %c:", 120+j);
                for (i=0;i<3;i++) {
                        printf(" % 1.7f", L[i][j]);
                }
                printf("\n");
        }
	printf("+--A1 2 3-----------------------------------------+\n");
	for (j=0;j<3;j++) {
		printf(" %c:", 120+j);
		for (i=0;i<3;i++) {
			printf(" % 1.7f", A[i][j]);
		}
		printf("\n");
	}
        printf("+--D1 2 3-----------------------------------------+\n");
        for (j=0;j<3;j++) {
		printf(" %c:", 120+j);
                for (i=0;i<3;i++) {
                        printf(" % 1.7f", D[i][j]);
                }
                printf("\n");
        }

	/* Calculate R */

	printf("\n*** Assumimg heliocentric orbit, geocentric pos'n. ***\n");
	printf("*** R will be read from the datafile (<< Almanac). ***\n\n");

	/* Calculate R_cos_psi = R_dot_L for t2 */

	R_cos_psi = dot(R[1],L[1]);

	/* Calculate R_dot_R for t2 */

	R_dot_R = dot(R[1],R[1]);

 	printf("+--R_cos_psi, R_dot_R-----------------------------+\n");
	printf("    % 1.7f, % 1.7f\n", R_cos_psi, R_dot_R);

	/* Calculate L_dot */
	/* p.384 */

	denom = -1*tau1*tau3*tau13;

 	printf("+-- -tau1.tau3.tau13------------------------------+\n");
	printf("    % 1.6f  % 1.6f  % 1.6f\n", tau1, tau3, tau13);
	printf("    % 1.7f\n", denom);

	G[0] = pow(tau3,2)/denom;
	G[2] = pow(tau1,2)/denom;
	G[1] = G[0] - G[2];
	L_dot[0] = -1*G[0]*L[0][0] + G[1]*L[1][0] + G[2]*L[2][0];
	L_dot[1] = -1*G[0]*L[0][1] + G[1]*L[1][1] + G[2]*L[2][1];
	L_dot[2] = -1*G[0]*L[0][2] + G[1]*L[1][2] + G[2]*L[2][2]; 

	R_dot[0] = -1*G[0]*R[0][0] + G[1]*R[1][0] + G[2]*R[2][0];
	R_dot[1] = -1*G[0]*R[0][1] + G[1]*R[1][1] + G[2]*R[2][1];
	R_dot[2] = -1*G[0]*R[0][2] + G[1]*R[1][2] + G[2]*R[2][2];

	/* Calculate L_dot_squared */

	L_dot_squared = dot(L_dot, L_dot);

	/* Calculate L_ddot, R_ddot */

	J[0] =  (2*tau3)/denom;
	J[2] = (-2*tau1)/denom;
	J[1] = J[0] + J[2];
	L_ddot[0] = J[0]*L[0][0] - J[1]*L[1][0] + J[2]*L[2][0];
	L_ddot[1] = J[0]*L[0][1] - J[1]*L[1][1] + J[2]*L[2][1];
	L_ddot[2] = J[0]*L[0][2] - J[1]*L[1][2] + J[2]*L[2][2];
	R_ddot[0] = J[0]*R[0][0] - J[1]*R[1][0] + J[2]*R[2][0];
	R_ddot[1] = J[0]*R[0][1] - J[1]*R[1][1] + J[2]*R[2][1];
	R_ddot[2] = J[0]*R[0][2] - J[1]*R[1][2] + J[2]*R[2][2];

	printf("+--G, J-------------------------------------------+\n");
	printf(" 1: % 1.6f,  % 2.6f\n", G[0], J[0]);
	printf(" 2: % 1.6f,  % 2.6f\n", G[1], J[1]);
	printf(" 3: % 1.6f,  % 2.6f\n", G[2], J[2]);

	printf("+--2L_dot, L_ddot, L_dot_squared------------------+\n");
	printf(" x: % 1.6f,  % 1.6f\n", 2*L_dot[0], L_ddot[0]);
	printf(" y: % 1.6f,  % 1.6f\n", 2*L_dot[1], L_ddot[1]);
	printf(" z: % 1.6f,  % 1.6f;  L_dot_squared: % 1.6f\n", \
	                 2*L_dot[2], L_ddot[2], L_dot_squared);

	/* Calculate E, A, B, A prime, B prime */
	/* p.378 */

	E_      =    determinant(L[1], L_dot, L_ddot);
	A_prime =    determinant(L[1], L_dot, R_ddot);
	B_prime = -1*determinant(L[1], L_dot, R[1]);
	A_ = A_prime/E_;
	B_ = B_prime/E_;

	printf("\n =[12B1]==========\n");
	printf("  E : % 1.7f\n  A': % 1.7f\n  B': % 1.7f\n",E_,A_prime,B_prime);
	printf(" =================\n");

	/* Estimate rho */
	/* p.386 */

	r_sq = estimrho(rho+1, R_cos_psi, R_dot_R, A_prime, B_prime, E_);

	printf(" rho[1]: % 1.7f\n",rho[1]);

	r_cu = pow(r_sq, -1.5);
	printf(" r^-3: % 1.7f", r_cu);

	/* Calculate approximate rho_dot [12E1] */
	/* pp. 391-393 */

	/* printf("\n\n"); */

	for (i=0;i<3;i++) {
		gin1[i][0] =   L[1][i];
		gin1[i][1] = 2*L_dot[i];
		gin1[i][2] =   L_ddot[i];
		gin1[i][3] =   R_ddot[i]+R[1][i]*r_cu;

		/* printf("\t% 1.8f % 1.8f % 1.8f % 1.8f\n", gin1[i][0],gin1[i][1],gin1[i][2],gin1[i][3]); */
	}

	gauss(&gin1[0][0], &rho_dot[0], 3);	/* :) */	/* Gaussian elimination using algorithm */
								/* from book by P. Wroblewski (Helios). */

	printf("\t rho_dot: % 1.7f\n", rho_dot[1]);

	/* Calculate new r, r_dot (first approx'n.) */
	/* p. 393 */

	printf("   \t\t .\n  r\t\t r\n --------------------------\n");
	for (i=0;i<3;i++) {
		r[1][i]     = rho[1]*L[1][i] - R[1][i];		
		r_dot[1][i] = rho_dot[1]*L[1][i] + rho[1]*L_dot[i] - R_dot[i];
		printf(" % 1.7f\t% 1.7f\n",r[1][i],r_dot[1][i]);
	}

	/* Observed time correction */
	/* p393 */

	tau1 = tau1*(1-k*alpha*rho_dot[1]);
	tau3 = tau3*(1-k*alpha*rho_dot[1]);

	/* f and g time series for r1, r3 */
	/* pp. 204-206 and 394 */

	g[0] = tau1;
	g[1] = tau3;

	for (i=0;i<3;i++) {
		ro        += r[1][i]*r[1][i];
		ro_dot    += r[1][i]*r_dot[1][i];
		so_dot_sq += r_dot[1][i]*r_dot[1][i];
	}
	ro     = sqrt(ro);
	ro_dot = ro_dot/ro;
	u      = fabs(so_dot_sq - 2/ro);

	fcoeff[0] = -0.5*pow(ro,-3);		/* f2, pp. 205-206 */
	gcoeff[0] = (1/3)*fcoeff[0];		/* g3 */
	fcoeff[1] = -0.5*ro_dot*pow(ro,-4);	/* f3 */
	gcoeff[1] = 0.5*fcoeff[1];		/* g4 */
	fcoeff[2] = (1/24)*pow(ro,-6)*(4 + (3*u*ro - 15*ro*ro_dot*ro_dot)); /* f4 */
	gcoeff[2] = 0.6*(fcoeff[2] - gcoeff[0]*gcoeff[0]); /* g5 */
	fcoeff[3] = -0.125*ro_dot*pow(ro,-7)*(4 + (3*u*ro - 7*ro*ro_dot*ro_dot)); /* f5 */
	gcoeff[3] = (2/3)*fcoeff[3] - gcoeff[0]*gcoeff[4]; /* g6 */
	
	for (i=2;i<6;i++) {
		f[0] += pow(tau1,i)*fcoeff[i-2];
		f[1] += pow(tau3,i)*fcoeff[i-2];
		g[0] += pow(tau1,i+1)*gcoeff[i-2];
		g[1] += pow(tau3,i+1)*gcoeff[i-2];
	}

	printf("\n           (1)\t\t  (3)\n ----------------------------------\n");
	printf("    tau: % 1.7f\t% 1.7f\n",tau1,tau3);
	printf("      f: % 1.7f\t% 1.7f\n      g: % 1.7f\t% 1.7f\n\n", f[0],f[1],g[0],g[1]);

	for (i=0;i<3;i++) {
		r[0][i] = f[0]*r[1][i] + g[0]*r_dot[1][i];
		r[2][i] = f[1]*r[1][i] + g[1]*r_dot[1][i];
		printf("    r %c: % 1.7f\t% 1.7f\n", rstr[i],r[0][i],r[2][i]);
	}

	printf("\n");
	for (i=0;i<3;i++) {
		rho_[0][i] = r[0][i] + R[0][i];
		rho_[2][i] = r[2][i] + R[2][i];
		printf("   rho %c: % 1.7f\t% 1.7f\n", rstr[i],rho_[0][i],rho_[2][i]);
	}

	/* Calculate residuals [12E3] */
	/* p. 394 */

	rho[0] = L[0][0]*rho_[0][0] + L[0][1]*rho_[0][1] + L[0][2]*rho_[0][2];
	rho[2] = L[2][0]*rho_[2][0] + L[2][1]*rho_[2][1] + L[2][2]*rho_[2][2];
        printf("\n  rho : % 1.7f\t% 1.7f\n",rho[0],rho[2]);

	cosd_dela[0] = -1*(A[0][0]*rho_[0][0] + A[0][1]*rho_[0][1])/rho[0];
	cosd_dela[2] = -1*(A[2][0]*rho_[2][0] + A[2][1]*rho_[2][1])/rho[2];
        printf("\n cddela: % 1.7f\t% 1.7f\n",cosd_dela[0],cosd_dela[2]);

	deld[0] = -1*(D[0][0]*rho_[0][0] + D[0][1]*rho_[0][1] + D[0][2]*rho_[0][2])/rho[0];
        deld[2] = -1*(D[2][0]*rho_[2][0] + D[2][1]*rho_[2][1] + D[2][2]*rho_[2][2])/rho[2];
        printf("  deld : % 1.7f\t% 1.7f\n",deld[0],deld[2]);

	do {	/* repeat corrections until residuals (above) are small enough */
	printf("\n ==== Correcting ... ======================(%d)== \n",ncorr+1);

	/* Calculate Leuschner differential correction */
	/* p. 396 */

	for(i=0;i<4;i++) {
		for(j=0;j<5;j++) {
			gin2[i][j] = 0;
		}
	}
	for(i=0;i<3;i++) {
		for(j=0;j<2;j++) {
			f_[j][i] = (f[j]*L[1][i] + aj[j]*r[1][i] + bj[j]*r_dot[1][i]);
		}
	}

	for(j=0;j<2;j++) {
		for (i=0;i<3;i++) {
			gin2[j][0]  +=A[j*2][i]*f_[j][i];
			gin2[j+2][0]+=D[j*2][i]*f_[j][i];
		}
		gin2[j][1]   = g[j]*A[j*2][0];
		gin2[j][2]   = g[j]*A[j*2][1];

		gin2[j+2][1] = g[j]*D[j*2][0];
		gin2[j+2][2] = g[j]*D[j*2][1];
		gin2[j+2][3] = g[j]*D[j*2][2];

		gin2[j][4]   = rho[j*2]*cosd_dela[j*2];
		gin2[j+2][4] = rho[j*2]*deld[j*2];				
	}
	

/*	for(i=0;i<4;i++) {
		for(j=0;j<5;j++) {
			printf(" % 1.7f",gin2[i][j]);
		}
		printf("\n");
	}
	printf("\n");
*/
	gauss(&gin2[0][0], &delrd[0], 4);

	printf("\n");
	for(i=0;i<3;i++) {
		delr[i] = L[1][i]*delrd[0];
		printf(" %1.7f delta_r %c: % 1.7f\n",delrd[i+1],rstr[i],delr[i]);
	}

	printf("   \t\t .\n     r(t2)\t r(t2)\t\t  rho(t2)\n");
	printf(" -------------------------------------------\n");
	for (i=0;i<3;i++) {
		r[1][i]     += delr[i];		
		r_dot[1][i] += delrd[i+1];
		delrdo[i]    = delrd[i+1];
		printf(" %c: % 1.7f\t% 1.7f\t% 1.7f\n",rstr[i],r[1][i],r_dot[1][i],rho[i]);
		
	}
	printf("\n");

	rho[0]+=dot(L[0],f_[0])*delrd[0] + g[0]*(dot(L[0],delrdo)); 
	rho[1]+=delrd[0];
	rho[2]+=dot(L[2],f_[1])*delrd[0] + g[1]*(dot(L[2],delrdo)); 
	
	for (i=0;i<3;i++) {
		printf(" rho t%d: % 1.7f\n",i,rho[i]);
	}

	/* 12E6 */

	/* Kepler */

	printf("\n");

	ro        = sqrt(dot(r[1],r[1]));
	so_dot_sq = dot(r_dot[1],r_dot[1]);
	u         = (200/ro - 100*so_dot_sq)*0.01;		/* = 1/a */
	u_ 	  = fabs(u);

	printf(" a: % 1.7f\tecosE: % 1.7f\tesinE: % 1.7f\tn: % 1.7f \n\n",1/u,1-ro*u,dot(r[1],r_dot[1])*sqrt(u_),sqrt(u_*u_*u_));

	E[0] = kepler(so_dot_sq*ro-1,dot(r[1],r_dot[1])*sqrt(u_),sqrt(u_*u_*u_)*tau1);
	E[1] = kepler(so_dot_sq*ro-1,dot(r[1],r_dot[1])*sqrt(u_),sqrt(u_*u_*u_)*tau3);

	printf("\n      E: % 1.7f\t% 1.7f\n",E[0],E[1]);	

	g[0] = (100*tau1 - pow(u,-1.5)*(100*E[0] - 100*sin(E[0])))*0.01;
	g[1] = (100*tau3 - pow(u,-1.5)*(100*E[1] - 100*sin(E[1])))*0.01;

	for(i=0;i<2;i++) {
		f[i] = (100 - 100*(1/(u*ro))*(1-cos(E[i])))*0.01;
	}
	printf("\n");
	for (i=0;i<3;i++) {
		rho_[0][i] = f[0]*r[1][i] + g[0]*r_dot[1][i] + R[0][i];
		rho_[2][i] = f[1]*r[1][i] + g[1]*r_dot[1][i] + R[2][i];
		printf("  rho %c: % 1.7f\t% 1.7f\n", rstr[i],rho_[0][i],rho_[2][i]);
	}
/*	rho[0] = L[0][0]*rho_[0][0] + L[0][1]*rho_[0][1] + L[0][2]*rho_[0][2]; */
	rho[0] = dot(L[0], rho_[0]);

/*	rho[2] = L[2][0]*rho_[2][0] + L[2][1]*rho_[2][1] + L[2][2]*rho_[2][2]; */
	rho[2] = dot(L[2], rho_[2]);

        printf("\n   rho : % 1.7f\t% 1.7f\n",rho[0],rho[2]);

	cosd_dela[0] = -1*(A[0][0]*rho_[0][0] + A[0][1]*rho_[0][1])/rho[0];
	cosd_dela[2] = -1*(A[2][0]*rho_[2][0] + A[2][1]*rho_[2][1])/rho[2];
        printf("\n cddela: % 1.7f\t% 1.7f\n",cosd_dela[0],cosd_dela[2]);

/*	deld[0] = -1*(D[0][0]*rho_[0][0] + D[0][1]*rho_[0][1] + D[0][2]*rho_[0][2])/rho[0]; */
	deld[0] = -1*dot(D[0], rho_[0])/rho[0];

/*      deld[2] = -1*(D[2][0]*rho_[2][0] + D[2][1]*rho_[2][1] + D[2][2]*rho_[2][2])/rho[2]; */
	deld[2] = -1*dot(D[2], rho_[2])/rho[2];

        printf("  deld : % 1.7f\t% 1.7f\n",deld[0],deld[2]);

	ncorr++;

	} while ((fabs(cosd_dela[0])>1e-5)||(fabs(cosd_dela[2])>1e-5)||(fabs(deld[0])>1e-5)||(fabs(deld[2])>1e-5));

	/* That's it! Calculate elements [12E8] */
	/* pp. 399-401 */

	printf("\n   Elements (after %d Leuschner corrections) ",ncorr);
	printf("\n+-----------------------------------------------+\n");

	Do      = dot(r[1],r_dot[1]);
	ro_dot  = Do/ro;

	E[0]    = 1 - ro*u;	/* ecosEo */
	E[1]    = Do*sqrt(u_);	/* esinEo */ 

	e       = sqrt(E[0]*E[0] + E[1]*E[1]);
	M       = (180/M_PI)*(((fabs(E[1])/E[1])*acos(E[0]/e) - E[1]));
	epsilon = rad(epsilon);
	coseps  = cos(epsilon);
	sineps  = sin(epsilon);

	/* printf("\t% 1.7f\t% 1.7f\n",E[0]/e,E[1]/e); */

	for (i=0;i<3;i++) {
		U_[i] = r[1][i]/ro;
		V_[i] = ro*r_dot[1][i] - Do*U_[i];
	}
	p = dot(V_,V_);
	for (i=0;i<3;i++) {
		V_[i] = V_[i]/sqrt(p);
		/* printf("\tU_[%d]: % 1.7f\tV_[%d]: % 1.7f\n", i,U_[i],i,V_[i]); */
	}
	Uzeo  = U_[2]*coseps-U_[1]*sineps;
	Vzeo  = V_[2]*coseps-V_[1]*sineps;

	incl  = fabs(asin(sqrt(pow(Uzeo,2)+pow(Vzeo,2))));

	dummy = (1/e)*ro_dot*sqrt(p);
	vo    = fabs(acos((1/e)*((p/ro) - 1)));
	vo    = (fabs(dummy)/dummy)*vo;

	dummy = (Uzeo/sin(incl));
	uo    = (fabs(dummy)/dummy)*fabs(acos(Vzeo/sin(incl)));

	dummy = (U_[1]*coseps + U_[2]*sineps - V_[0])/(1 + cos(incl));
	lo    = fabs(acos((U_[0] + V_[1]*coseps + V_[2]*sineps)/(1 + cos(incl))));
	lo    = (fabs(dummy)/dummy)*lo;

	printf("|  parameter     p: % 1.7f\t[% 1.7f]\t|\n",p,(1/u)*(1 - e*e));
	printf("|                a: % 1.7f\t\t\t|\n",1/u);
	printf("|  eccentricity  e: % 1.8f\t\t\t|\n",e);	
	printf("|  mean anomaly  M: % 1.7f deg\t\t|\n",M);
	printf("|  true anomaly  v: % 3.7f deg\t\t|\n",vo*180.0/M_PI);
	printf("|                u: % 3.7f deg\t\t|\n",uo*180.0/M_PI);
	printf("|                l: % 3.7f deg\t\t|\n",lo*180.0/M_PI);
	printf("|  inclination   i: % 3.7f deg\t\t|\n",incl*180.0/M_PI);
	printf("|  Omega         O: % 3.7f deg\t\t|\n",360 + (lo - uo)*180.0/M_PI);
	printf("|  omega         w: % 3.7f deg\t\t|\n",(uo - vo)*180.0/M_PI); 
	printf("|  n(t)          n: % 1.9f deg\t\t|\n",(180/M_PI)*pow(u,1.5)*k);
	printf("+-----------------------------------------------+\n\n");


	/* --------------------------------------------------- */

	exit(0);
}

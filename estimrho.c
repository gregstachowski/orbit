#include <math.h>
#include <stdio.h>

double calc_f_rho(double A, double B, double r_sq, double rho) {
	return(A - B*pow(r_sq,-1.5) - rho);
}	

double estimrho(	double *rhop,
			double R_cos_psi,
			double R_dot_R,
			double A_prime,
			double B_prime,
			double E_	) {

	double A_ = A_prime/E_;
	double B_ = B_prime/E_;
	double f_rho, f_prime, delta_rho;
	double rho1, rho2, r_sq1, r_sq2;

	printf("\n");
	rho1 = 3.0;
	rho2 = 10.0;
	
	do {
		/* printf(" % 0.6f  % 0.6f",rho1, rho2); */
		r_sq1 = (rho1*rho1 - 2*R_cos_psi*rho1 + R_dot_R);
		r_sq2 = (rho2*rho2 - 2*R_cos_psi*rho2 + R_dot_R);
		f_rho = calc_f_rho(A_,B_,r_sq1,rho1);
		f_prime = (f_rho-calc_f_rho(A_,B_,r_sq2,rho2))/(rho1-rho2);
		delta_rho = -1*f_rho/f_prime;
		/* printf("\t% 1.9f\t% 1.9f\t% 1.7f\n",delta_rho, f_rho, f_prime); */
		rho1 += delta_rho*0.25;
	} while (fabs(f_rho) > 1e-9); 	
	printf("\n");

	*rhop = rho1;
	return(r_sq1);
}

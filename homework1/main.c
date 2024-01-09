#include <stdio.h>
#include "nr.h"

void get_eps(float *eps){
	*eps = 0.f;
	float epsilon_float = 1.f;

	while (1) {
		*eps = epsilon_float;
		epsilon_float /= 2.f;

		if (1.f + epsilon_float == 1.f) {
			break;
		}
	}
}

void get_eps_double(double *eps){
	*eps = 0.0;
	double epsilon_double = 1.0;

	while (1) {
		*eps = epsilon_double;
		epsilon_double /= 2.0;

		if (1.0 + epsilon_double == 1.0) {
			break;
		}
	}
}

int main(){
	int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	float eps, epsneg, xmin, xmax;
	double eps1, epsneg1, xmin1, xmax1;

	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp,
			&eps, &epsneg, &xmin, &xmax);
	printf("Machine Accuracy (machar): \t%0.20f\n", eps);

	get_eps(&eps);
	printf("Machine Accuracy (get_eps): \t%0.20f\n", eps);

	machar1(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp,
			&eps1, &epsneg1, &xmin1, &xmax1);
	printf("Machine Accuracy (machar1): \t%0.20f\n", eps1);

	get_eps_double(&eps1);
	printf("Machine Accuracy (get_eps_d): \t%0.20f\n", eps1);
	
	return 0;
}

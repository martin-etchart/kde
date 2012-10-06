#include <stdio.h>
#include <math.h>
#include "fftw3.h"

#define N 3

int main(int argc, char** argv)
{
	int n=N;

	double in[N]={2.0,1.0,0.0};
	double out[N], in2[N];

	fftw_plan idct =	fftw_plan_r2r_1d(3, in, out, FFTW_REDFT01,
			FFTW_MEASURE);
	fftw_plan dct =	fftw_plan_r2r_1d(3, out, in2, FFTW_REDFT10,
			FFTW_MEASURE);


	for(int i = 0; i < n; i++)
		in[i]*=1.0/sqrt(2*n);

	in[0]*=sqrt(2);


	fftw_execute(idct);
	fftw_execute(dct);

	for(int i = 0; i < n; i++)
		in2[i]*=sqrt(2*n);

	printf("in: ");
	for(int i = 0; i < n; i++)
	{
		printf(" %g", in[i]);
	}
	printf("\n");

	printf("out: ");
	for(int i = 0; i < n; i++)
	{
		printf(" %g", out[i]);
	}
	printf("\n");

	printf("in2: ");
	for(int i = 0; i < n; i++)
	{
		printf(" %g", in2[i]);
	}
	printf("\n");


	fftw_destroy_plan(idct);
	fftw_destroy_plan(dct);

	return 0;
}

#include <stdio.h>
#include <math.h>
#include "fftw3.h"

#define N 3

void print_v(double* v, int n, char* title)
{
	printf("%s:", title);

	for(int i = 0; i < n; i++)
	{
		printf(" %g", v[i]);
	}
	printf("\n");
}

int main(int argc, char** argv)
{
	int n=N;

	double in[N]={2.0,1.0,0.0};
	double out[N],in_n[N], in2[N],in2_n[N];

	fftw_plan idct =	fftw_plan_r2r_1d(N, in_n, out, FFTW_REDFT01,
			FFTW_MEASURE);
	fftw_plan dct =	fftw_plan_r2r_1d(N, out, in2, FFTW_REDFT10,
			FFTW_MEASURE);


	for(int i = 0; i < n; i++)
		in_n[i]=in[i]/sqrt(2*n);

	in_n[0]*=sqrt(2);


	fftw_execute(idct);
	fftw_execute(dct);

	for(int i = 0; i < n; i++)
		in2_n[i]=in2[i]/sqrt(2*n);

	in2_n[0]/=sqrt(2);

	print_v(in,n,"in");
	print_v(in_n,n,"in_n");
	print_v(out,n,"out");
	print_v(in2,n,"in2");
	print_v(in2_n,n,"in2_n");

	fftw_destroy_plan(idct);
	fftw_destroy_plan(dct);

	return 0;
}

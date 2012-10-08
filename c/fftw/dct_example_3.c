#include <stdio.h>
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
	double out[N], in2[N],in2_n[N];

	fftw_plan idct =	fftw_plan_r2r_1d(3, in, out, FFTW_REDFT01,
			FFTW_MEASURE);
	fftw_plan dct =	fftw_plan_r2r_1d(3, out, in2, FFTW_REDFT10,
			FFTW_MEASURE);


	fftw_execute(idct);
	fftw_execute(dct);


	for(int i = 0; i < 3; i++)
		in2_n[i]=in2[i]/(2*n);


	print_v(in,N,"in");
	print_v(out,N,"out");
	print_v(in2,N,"in2");
	print_v(in2_n,N,"in2_n");



	fftw_destroy_plan(idct);
	fftw_destroy_plan(dct);

	return 0;
}

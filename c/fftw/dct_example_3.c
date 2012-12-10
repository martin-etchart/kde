#include <stdio.h>
#include "fftw3.h"



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
#define N 5

	double in[N]={2.00,0,-0.93, 0, 0.08};
	//double in[N]={2.0,1.0,0.0};
	//double in[N]={1.41 ,0,-0.82, 0, 0.16,0};
	//double in[N]={0.0545,0.2442,0.4026,0.2442,0.0545};
	double out[N];
	//double out[N]={0.017559500000000 ,  0.129748000000000  , 0.352692000000000   ,0.352692000000000  , 0.129748000000000  , 0.017559500000000};
	//double out[N]={0.017559513479670,   0.129748230171210,   0.352692256349121,   0.352692256349121,   0.129748230171210,   0.017559513479670};
	//double out[N]={0.0545,0.2442,0.4026,0.2442,0.0545};
	double in2[N],in2_n[N];
	int n=N;

	fftw_plan idct =	fftw_plan_r2r_1d(n, in, out, FFTW_REDFT01, FFTW_MEASURE);
	fftw_plan dct =	fftw_plan_r2r_1d(n, out, in2, FFTW_REDFT10, FFTW_MEASURE);


	fftw_execute(idct);
	fftw_execute(dct);


	for(int i = 0; i < n; i++)
		in2_n[i]=in2[i]/(2*n);


	print_v(in,N,"in");
	print_v(out,N,"out");
	print_v(in2,N,"in2");
	print_v(in2_n,N,"in2_n");



	fftw_destroy_plan(idct);
	fftw_destroy_plan(dct);

	return 0;
}

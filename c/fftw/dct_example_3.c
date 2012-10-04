#include <stdio.h>
#include "fftw3.h"

#define N 3

int main(int argc, char** argv)
{

int n=N;

double in[3]={2.0,1.0,0.0};
double out[3], in2[3];

fftw_plan idct =	fftw_plan_r2r_1d(3, in, out, FFTW_REDFT01,
FFTW_MEASURE);
fftw_plan dct =	fftw_plan_r2r_1d(3, out, in2, FFTW_REDFT10,
FFTW_MEASURE);


fftw_execute(idct);
fftw_execute(dct);


for(int i = 0; i < 3; i++)
 in2[i]/=2*n;


for(int i = 0; i < 3; i++)
{
printf(" %g", out[i]);
}
printf("\n");


for(int i = 0; i < 3; i++)
{
printf(" %g", in2[i]);
}
printf("\n");


fftw_destroy_plan(idct);
fftw_destroy_plan(dct);

return 0;
}

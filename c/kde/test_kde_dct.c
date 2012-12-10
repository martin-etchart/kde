#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#include <complex.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "kde_util.h"

#include <fftw3.h>

int verbose = -1;

void print_v(double* v, int n, char* title)
{
	printf("%s:", title);

	for(int i = 0; i < n; i++)
	{
		printf(" %g", v[i]);
	}
	printf("\n");
}


int main( int argc, char** argv )
{
	XML_IN;
	int length=0;
	double *data=NULL;
	const char * full_fname = "../../../matlab/test_data.txt";

	file_read_into_array_doubles(full_fname, &data, &length);


	double *out=malloc(length*sizeof(*out));
	kde_dct_fftw(data, length, out);
	dct_fftw(data, length, out);

	double *data_rec=malloc(length*sizeof(*data_rec));
	kde_idct_fftw(data, length, data_rec);
	idct_fftw(data, length, data_rec);

	if  (verbose==1 || verbose==-1)
	{
		print_vec(data,"data",0,length);
		print_vec(out,"out",0,length);
		print_vec(data,"data_rec",0,length);
	}

	free(out);
	free(data);
	free(data_rec);
	XML_OUT;
	return 0;
}

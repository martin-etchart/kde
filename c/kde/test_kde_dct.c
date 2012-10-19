#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc/malloc.h>
#include <complex.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "kde_util.h"

int verbose = -1;



int main( int argc, char** argv )
{
	XML_IN;
	int length=0;
	double *data=NULL;
	const char * full_fname = "../../../matlab/test_data.txt";

	file_read_into_array_doubles(full_fname, &data, &length);


	double *out=malloc(length*sizeof(*out));
	kde_dct_fftw(data, length, out);

	if  (verbose==1 || verbose==-1)
	{
		print_vec(data,"data",0,length);
		print_vec(out,"out",0,length);
	}

	free(out);
	free(data);
	XML_OUT;
	return 0;
}

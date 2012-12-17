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
#include "kde.h"

#include <fftw3.h>


int main( int argc, char** argv )
{
	XML_IN;
	int length=0;
	int length_x=0;
	double *data=NULL;
	double *x=NULL;
	const char * full_fname = "../../../data/test/test_extrema.txt";
	const char * full_fname_x = "../../../data/test/test_extrema_x.txt";

	file_read_into_array_doubles(full_fname, &data, &length);
	file_read_into_array_doubles(full_fname_x, &x, &length_x);

	double delta=1e-3;

	int l_min,l_max;
	double* min_x;
	double* max_x;
	peakdet( length, x, data, delta,&l_min,&min_x,&l_max,&max_x);

	if  (verbose==1 || verbose==-1)
	{
		print_vec(x,"x",0,length);
		print_vec(data,"data",0,length);
		print_vec(min_x,"min_x",0,l_min);
		print_vec(max_x,"max_x",0,l_max);
	}

	free(data);

	XML_OUT;
	return 0;
}

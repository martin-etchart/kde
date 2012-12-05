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



int main( int argc, char** argv )
{
XML_IN;
	int length=0;
	double *data=NULL;
	const char * full_fname = "../../../matlab/data.txt";
	
	file_read_into_array_doubles(full_fname, &data, &length);
	
	if  (verbose==1 || verbose==-1)
	{
	print_vec(data,"data",0,length);
	}
	
	
	if (verbose==1 || verbose==-1)
	{
		//printf("---DATA---\n"); for (int i=0; i<300; i++) printf("%f\n",data[i]);
	}
	free(data);
XML_OUT;
	return 0;
}

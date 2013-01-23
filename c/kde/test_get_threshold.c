#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#include "kde_util.h"
#include "roots.h"
#include "kde.h"



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

	double maximum, minimum;
	double bw=-1;
	double *density=NULL;
	double *x=NULL;
	find_max_min_array_doubles(data,length,&maximum,&minimum);
	int n=128;
	kde(data,length,n,minimum-5,maximum+5, &density, &x, &bw);

	//compute maxima
	double delta=1e-3;
	int l_min,l_max;
	double* min_x;
	double* max_x;
	peakdet( n, x, density, delta, &l_min,&min_x,&l_max,&max_x);

	if  (verbose==1 || verbose==-1)
	{
		print_vec(x,"x",0,n);
		print_vec(density,"density",0,n);
		print_vec(data,"data",0,length);
		print_vec(min_x,"min_x",0,l_min);
		print_vec(max_x,"max_x",0,l_max);
	}



	if(!density)
		free(density);
	free(data);
	XML_OUT;
	return 0;
}

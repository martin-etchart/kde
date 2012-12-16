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
	find_max_min_array_doubles(data,length,&maximum,&minimum);
	//kde(data,length,pow(2,14),minimum-5,maximum+5);
	kde(data,length,128,minimum-5,maximum+5);



	free(data);
	XML_OUT;
	return 0;
}

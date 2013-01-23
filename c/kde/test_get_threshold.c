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
	//const char * full_fname = "../../../matlab/data.txt";
	const char * full_fname = "./data_in.txt";

	file_read_into_array_doubles(full_fname, &data, &length);

	bones_get_threshold(data,length);


	free(data);
	XML_OUT;
	return 0;
}

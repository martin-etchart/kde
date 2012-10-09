#include "kde_util.h"
#include <stdio.h>

int print_vec(double *v, char* title, int start, int end)
{
	XML_IN_T(title);
	//printf("%s: ",title);
	for(int i=start;i<end;i++)
		printf(" %g",v[i]);
	printf("\n");
	XML_OUT;
}

int file_read_into_array_doubles(const char *filename , double **out_data, int *length)
{
    FILE *in_file;
    in_file = fopen(filename, "r");
	 double *data=NULL;

    if (in_file == NULL)
    {
        return -1;
    }
    else
    {
		 fscanf(in_file, "%d", length);
		 printf("length: %d\n", *length);
		 data =(double*)malloc((*length)*sizeof(*data));

        for(int j=0; j<*length; j++)
        {
            fscanf(in_file, "%lg", data+j);
        }
        fclose(in_file);
		  *out_data=data;
    }
    return 0;
}


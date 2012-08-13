#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_histogram.h>

int verbose = -1;

int file_read_into_array_doubles(const char *filename , double *data, int *length)
{
    FILE *in_file;
    in_file = fopen(filename, "r");

    if (in_file == NULL)
    {
        return -1;
    }
    else
    {
        for(int j=0; j<*length; j++)
        {
            fscanf(in_file, "%lf", &data[j]);
        }
        fclose(in_file);
    }
    return 0;
}

int find_max_min_array_doubles(double *a, int length, double *max, double *min)
{
	// Find maximum and minimum element in an arary.
	*max = a[0];
	*min = a[0];

	for (int i = 0; i < length; i++)
	{
		if (a[i] > *max)
			*max = a[i];
		else if (a[i] < *min)
			*min = a[i];
	}
	if  (verbose==2 || verbose==-1)
	{
		printf ("Maximum element in an array : %f\n", *max);
		printf ("Minimum element in an array : %f\n", *min);
	}

 	return 0;
}

int compare_doubles (const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;
	
	return (*da > *db) - (*da < *db);
}

double log2(double x) {	return (log(x) / log(2.0)); }

int histc(double *data, double *xmesh, int n , double *bins){

	printf("-1\n");

	gsl_histogram * h = gsl_histogram_alloc(n-1);

	printf("0\n");

	gsl_histogram_set_ranges (h, xmesh, n);

	printf("1\n");

	for (int i=0; i<300; i++)
		gsl_histogram_increment (h, data[i]);

	printf("2\n");

	for (int i=0;i<300;i++)
		bins[i] = gsl_histogram_get (h, i);

	printf("3\n");

	gsl_histogram_fprintf (stdout, h, "%g", "%g");
	/*...*/

	//gsl_histogram_free(h);

	return 0;
}

void kde(double *data, int length, int n ,double dataMIN, double dataMAX)
{
/*
 * function [bandwidth,density,xmesh,cdf]=kde(data,n,dataMIN,dataMAX)
 *  Reference: 
 * Kernel density estimation via diffusion
 * Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
 * Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
*/
	
	/* if n is not supplied switch to the default
	n = pow( 2 , 14 ); 
	*/
	n = pow( 2 , ceil(log2(n)) ); // round up n to the next power of 2;
	/*define the default  interval [MIN,MAX]
	double maximum,minimum;
	find_max_min_array_doubles(data,length,&maximum,&minimum);
	double Range=maximum-minimum;
    dataMIN=minimum-Range/10; dataMAX=maximum+Range/10;
    */
    
    /*set up the grid over which the density estimate is computed;*/
	double R=dataMAX-dataMIN; 
	double dx=R/(n-1);
	double xmesh[n];
	for (int i=0; i<n; i++ ) xmesh[i] = dataMIN+i*dx;

	//...double N=length(unique(data));
	//qsort(data, length, sizeof(double), compare_doubles);
	
	/*bin the data uniformly using the grid defined above;*/
	double initial_data[n];

	printf("a\n");
	histc(data, xmesh, n , initial_data);
	printf("b\n");
//	initial_data=histc(data,xmesh)/N;  
//	initial_data=initial_data/sum(initial_data);
//	a=dct1d(initial_data); // discrete cosine transform of initial data
	
	/*now compute the optimal bandwidth^2 using the referenced method*/
//	I=[1:n-1]'.^2; 
//	a2=(a(2:end)/2).^2;
	
	/*use  fzero to solve the equation t=zeta*gamma^[5](t)*/
/*	try
		t_star=fzero(@(t)fixed_point(t,N,I,a2),[0,.1])
	catch
		t_star=.28*N^(-2/5);
	end
*/	
	/*smooth the discrete cosine transform of initial data using t_star*/
//	a_t=a.*exp(-[0:n-1]'.^2*pi^2*t_star/2);
	
	
	
	if  (verbose==3 || verbose==-1)
	{
	for (int i=1700; i<2000 ; i++ )
		//printf("%f\n",xmesh[i]);
		printf("%f\n",initial_data[i]);
	}

	
}


int main( int argc, char** argv )
{

	int length = 300;
	double data[length];
	const char * full_fname = "/home/roho/workspace/juanc/matlab/matlab/kde/data.txt";
	
	file_read_into_array_doubles(full_fname, data, &length);
	
	if  (verbose==1 || verbose==-1)
	{
		//printf("---DATA---\n"); for (int i=0; i<300; i++) printf("%f\n",data[i]);
	}
	
	double maximum, minimum;
	find_max_min_array_doubles(data,length,&maximum,&minimum);
	kde(data,length,pow(2,14),minimum-5,maximum+5);
	
	if (verbose==1 || verbose==-1)
	{
		//printf("---DATA---\n"); for (int i=0; i<300; i++) printf("%f\n",data[i]);
	}

	return 0;
}

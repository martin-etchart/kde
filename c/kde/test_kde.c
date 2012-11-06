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


int file_read_into_array_doubles_l(const char *filename , double *data, int *length)
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

int histc(double *data, double *xmesh, int n , double *bins)
{

	printf("-1\n");

	gsl_histogram * h = gsl_histogram_alloc(n-1);

	printf("0\n");

	gsl_histogram_set_ranges (h, xmesh, n);
	
	double h_max = gsl_histogram_max(h);
	double h_min = gsl_histogram_min(h);

	printf("1\n");

	for (int i=0; i<300; i++)
		gsl_histogram_increment (h, data[i]);

	printf("2\n");

	for (int i=0;i<n-1;i++)
		bins[i] = gsl_histogram_get (h, i);

	printf("3\n");

	gsl_histogram_fprintf (stdout, h, "%g", "%g");
	/*...*/

	gsl_histogram_free(h);

	return 0;
}

int fft(double *data, int length, double complex *fft_data)
{

	gsl_fft_real_radix2_transform (data, 1, length);

	return 0;
}

int ifft(double *data, int length, double *ifft_data)
{


	gsl_fft_halfcomplex_radix2_inverse (data, 1, length);

	return 0;
}

int dct1d(double *data, int length, double *dct_data)
{
/*	% computes the discrete cosine transform of the column vector data
	[nrows,ncols]= size(data);
	% Compute weights to multiply DFT coefficients
	weight = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
	% Re-order the elements of the columns of x
	data = [ data(1:2:end,:); data(end:-2:2,:) ];
	% Multiply FFT by weights:
	data= real(weight.* fft(data));
*/
	/*Compute weights to multiply DFT coefficients*/
	double complex weight[length];
	weight[0] = 1;
	for (int i=1;i<length;i++)
		weight[i] = 2*(cexp(-I*i*M_PI/(2*length)));

	/*Re-order the elements of the columns of x*/
	double data_aux[length/2];
	for(int i=0, j=length-1 ; j>=2 ; i++, j-=2 )
		data_aux[i] = data[j];

	/*Multiply FFT by weights*/
	double complex fft_data[length];
	fft(data, length, fft_data);
	for (int i=0;i<length;i++)
		dct_data[i] = creal(weight[i]*fft_data[i]);

	return 0;
}

int idct1d(double *data, int length, double *dct_data)
{
/*	% computes the inverse discrete cosine transform
	[nrows,ncols]=size(data);
	% Compute weights
	weights = nrows*exp(i*(0:nrows-1)*pi/(2*nrows)).';
	% Compute x tilde using equation (5.93) in Jain
	data = real(ifft(weights.*data));
	% Re-order elements of each column according to equations (5.93) and
	% (5.94) in Jain
	out = zeros(nrows,1);
	out(1:2:nrows) = data(1:nrows/2);
	out(2:2:nrows) = data(nrows:-1:nrows/2+1);
	%   Reference:
	%      A. K. Jain, "Fundamentals of Digital Image
	%      Processing", pp. 150-153.
*/

	/*Compute weights*/
	double complex weight[length];
	for (int i=0;i<length;i++)
		weight[i] = length*(cexp(I*i*M_PI/(2*length)));

	/*Compute weighted data*/
	for (int i=0;i<length;i++)
		data[i] = data[i]*weight[i];

	/*get IFFT of weighted data*/
	double complex ifft_data[length];
	ifft(data, length, ifft_data);

	/*Compute x tilde using equation (5.93) in Jain*/
	for (int i=0;i<length;i++)
		dct_data[i] = creal(ifft_data[i]);

	/*Re-order elements of each column according to equations (5.93) and (5.94) in Jain*/
	//out = zeros(nrows,1);
	//out(1:2:nrows) = data(1:nrows/2);
	//out(2:2:nrows) = data(nrows:-1:nrows/2+1);

	return 0;
}

double fixed_point(double t, int N, double *It, double *a2, int n)
{
/*	function  out=fixed_point(t,N,I,a2)
	% this implements the function t-zeta*gamma^[l](t)
	*/
	//l=7;
	int l=7;
double sum_f=0;
for(int i=0;i<n;i++)
	sum_f+=pow(It[i],l)*a2[i]*exp(-It[i]*pow(MATH_PI,2.0)*t);
double f=2*pow(MATH_PI,2.0*l)*sum_f;

/*
f=2*pi^(2*l)*sum(I.^l.*a2.*exp(-I*pi^2*t));
	for s=l-1:-1:2
		K0=prod([1:2:2*s-1])/sqrt(2*pi);  const=(1+(1/2)^(s+1/2))/3;
		time=(2*const*K0/N/f)^(2/(3+2*s));
		f=2*pi^(2*s)*sum(I.^s.*a2.*exp(-I*pi^2*time));
	end
	out=t-(2*N*sqrt(pi)*f)^(-2/5);
	end
*/

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

	double N = length; //double N=length(unique(data)); //qsort(data, length, sizeof(double), compare_doubles);
	
	/*bin the data uniformly using the grid defined above;*/
	double initial_data[n];
	histc(data, xmesh, n , initial_data);
	double sum_initial_data = 0;
	for(int i=0;i<n;i++){
		initial_data[i]=initial_data[i]/N;
		sum_initial_data+=initial_data[i];
	}
	for(int i=0;i<n;i++)
		initial_data[i]=initial_data[i]/sum_initial_data;

	double a[n];
	kde_dct_fftw(initial_data,n,a); // discrete cosine transform of initial data
	
	/*now compute the optimal bandwidth^2 using the referenced method*/
	double It[n-1];
	for (int i=0;i<n-1;i++)
		It[i]=pow(i+1,2);
	double a2[n-1];
	for(int i=0;i<n-1;i++)
		a2[i] = pow(a[i+1]/2,2);

double t_star=.28*pow(N,-2.0/5.0);
	/*use  fzero to solve the equation t=zeta*gamma^[5](t)*/
/*	try
		t_star=fzero(@(t)fixed_point(t,N,I,a2),[0,.1])
	catch
		t_star=.28*N^(-2/5);
	end
*/
  
	/*smooth the discrete cosine transform of initial data using t_star*/
	double a_t[n];
	for(int i=0;i<n;i++)
		a_t[i]=a[i]*exp(-i^2*M_PI^2*t_star/2);
	
	
	
	if  (verbose==3 || verbose==-1)
	{
		print_vec(xmesh,"xmesh",0,128);
		print_vec(initial_data,"initial_data",0,128);
		print_vec(a,"a",0,128);
	}

	
}


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

	if (verbose==1 || verbose==-1)
	{
		//printf("---DATA---\n"); for (int i=0; i<300; i++) printf("%f\n",data[i]);
	}
	free(data);
	XML_OUT;
	return 0;
}

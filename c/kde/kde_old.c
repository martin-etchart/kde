#include <complex.h>

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



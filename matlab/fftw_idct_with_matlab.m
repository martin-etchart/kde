function idct_fftw=fftw_idct_with_matlab(x)

n=length(x);

x(1) = x(1) / sqrt(2);
idct_fftw= idct(x) /sqrt(2*n);
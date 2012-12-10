function idct_fftw=fftw_idct_with_matlab(x)

n=length(x);

xt=x;
xt(1) = xt(1)/sqrt(2);
idct_fftw= idct(xt)/sqrt(2*n) ;

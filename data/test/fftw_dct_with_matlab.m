function dct_fftw=fftw_dct_with_matlab(x)

n=length(x);

xt=x;
dct_fftw=dct(xt);
dct_fftw=2*dct_fftw*sqrt(n/2);
dct_fftw(1)=dct_fftw(1)*sqrt(2);

end
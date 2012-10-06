clear all
close all
home
%% make matlab give the same output as fftw
%input
x=[2 1 0]

%% regular matlab output
idct_matlab=idct(x)

%% regular dct
dct_matlab=dct(idct_matlab)

%% fftw dct
%normalized input
y = x;
y(1) = x(1) / sqrt(2);
idct_fftw= idct(y) *sqrt(2*length(x))

dct_fftw=dct(idct_fftw/sqrt(2*length(x)));
dct_fftw(1)=dct_fftw(1)*sqrt(2);
dct_fftw
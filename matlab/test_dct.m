clear all
close all
home
%% make matlab give the same output as fftw
%input
x=[2 1 0]
n=length(x);
%% regular matlab output
% idct
idct_matlab=idct(x)
% dct
dct_matlab=dct(idct_matlab)

%% fftw dct
%normalized input
y = x;
y(1) = x(1) / sqrt(2);
idct_fftw= idct(y) *sqrt(2*n)

xt=idct_fftw/sqrt(2*n)
dct_fftw=dct(xt);
dct_fftw(1)=dct_fftw(1)*sqrt(2);
dct_fftw

%% test direct DCT
in=[1.8618    1.1547    0.4476]
out_matlab=dct(in)
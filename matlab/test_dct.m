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
idct_fftw=fftw_idct_with_matlab(x)

dct_fftw=fftw_dct_with_matlab(idct_fftw)

%% test direct DCT
in=[1.8618    1.1547    0.4476]
out_matlab=dct(in)


%% test
%in=[0.0545    0.2442    0.4026    0.2442    0.0545];
in=[1.86181 1.1547 0.447594];
out_kde_dct=kde_dct1d(in)
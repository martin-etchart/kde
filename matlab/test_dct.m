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
in=[0.0545    0.2442    0.4026    0.2442    0.0545];
out_kde_dct=kde_dct1d(in)

%% test dct fftw using matlab
in=[0.0545    0.2442    0.4026    0.2442    0.0545];
out_kde_dct=kde_dct1d(in)

%% full test dct reproduced in matlab
in = fspecial('gaussian',[128 1],5);

out_fftw_matlab=fftw_dct_with_matlab(in);
norm(out_fftw_matlab-out_fftw')
figure(4)
plot(out_fftw_matlab,'-r')
hold on
plot(out_fftw,'-b')
hold off

%% full test idct reproduced in matlab

in_rec_fftw_matlab=fftw_idct_with_matlab(out_fftw);
norm(in_rec_fftw_matlab-in_rec_fftw)
figure(4)
plot(in_rec_fftw_matlab,'-r')
hold on
plot(in_rec_fftw,'-b')
hold off
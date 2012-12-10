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

%% fftw dct made with matlab
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

%% test dct fftw using matlab
%in=[0.0545    0.2442    0.4026    0.2442    0.0545];
in =[0.017559513479670   0.129748230171210   0.352692256349121   0.352692256349121   0.129748230171210   0.017559513479670];
out_fftw_dct=fftw_dct_with_matlab(in)

%% test idct fftw using matlab
%out=[2 1 0];
out=[ 2.00             0         -0.93             0          0.08];
in_rec_fftw_dct=fftw_idct_with_matlab(out)

%% full test dct reproduced in matlab
in = fspecial('gaussian',[128 1],10);
%in =[0.017559513479670   0.129748230171210   0.352692256349121   0.352692256349121   0.129748230171210   0.017559513479670]';
n=length(in);

if 1
	save('test_data','in')
	jfc_vector_save_simple(in,'test_data.txt');
end

%out_fftw=jfc_vector_read_simple('out_fftw.txt');
out_fftw_matlab=fftw_dct_with_matlab(in);
norm(out_fftw_matlab-out_fftw')
figure(41)
plot(out_fftw_matlab,'-*b')
hold on
plot(out_fftw,'-*r')
hold off

figure(42)
plot(in,'-b')
hold on
plot(in_fftw,'-r')
hold off
norm(in-in_fftw')

%% full test idct reproduced in matlab

in_rec_fftw_matlab=fftw_idct_with_matlab(out_fftw_matlab);
norm(in_rec_fftw_matlab-in_rec_fftw')
figure(4)
plot(in_rec_fftw_matlab,'-b')
hold on
plot(in_rec_fftw,'-r')
hold off
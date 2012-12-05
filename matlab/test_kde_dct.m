N=length(initial_data);

%%use kde_dct to reproduce matlab
in_data_n=initial_data;
in_data_n=in_data_n/sqrt(2*N);
y1=dct(initial_data);
y2=kde_dct1d(in_data_n);
y2(1)=y2(1)*sqrt(2);
figure(1)
plot(y1,'-r')
hold on
plot(y2,'-b')
hold off

dct_error.y=abs(y1-y2)./abs(y1);
dct_error.mean=mean(dct_error.y);
dct_error.l2=norm(dct_error.y);
dct_error.max=max(dct_error.y);

dct_error

figure(2)
plot(dct_error.y)

%% use matlab to reproduce kde_dct
% in_data_n=initial_data*sqrt(2*N);
% z1=dct(in_data_n);
% z1(1)=z1(1)/sqrt(2);
in = fspecial('gaussian',[1 64],2);
%in = initial_data;
z1=kde_dct_with_matlab(in);
z2=kde_dct1d(in);

figure(3)
plot(z1,'-r')
hold on
plot(z2,'-b')
hold off

dct_error.y=abs(z1-z2+eps)./abs(z1+eps);
dct_error.mean=mean(dct_error.y);

dct_error.max=max(dct_error.y);

dct_error

figure(4)
plot(dct_error.y)

%%

%%
%in=[1.8618    1.1547    0.4476];
%in=[0.0545    0.2442    0.4026    0.2442    0.0545]';
in = fspecial('gaussian',[128 1],5);

out_kde_dct=kde_dct1d(in);
out_matlab=kde_dct_with_matlab(in);
norm(out_kde_dct-out_matlab)
figure(4)
plot(out_kde_dct,'-r')
hold on
plot(out_matlab,'-b')
hold off

%% test kde vs fftw direct
in = fspecial('gaussian',[128 1],5);
n=length(in);

if 1
	save('test_data','in')
	jfc_vector_save_simple(in,'test_data.txt');
end
	
out_kde_dct=kde_dct1d(in);
in_n=in*2*n;
out_fftw=fftw_dct_with_matlab(in_n);
out_fftw(1)=out_fftw(1)/2;
norm(out_kde_dct-out_fftw)
figure(4)
plot(out_kde_dct,'-r')
hold on
plot(out_fftw,'-b')
hold off

%% test kde vs fftw inverse
in_kde_rec=kde_idct1d(out_kde_dct);
in_kde_rec_n=in_kde_rec/sum(in_kde_rec);
out_kde_dct_n=out_kde_dct;
out_kde_dct_n(1)=out_kde_dct_n(1)*2;
in_kde_rec_fftw=fftw_idct_with_matlab(out_kde_dct_n)/2;
in_kde_rec_fftw_n=in_kde_rec_fftw/sum(in_kde_rec_fftw);
norm(in_kde_rec_n-in)
norm(in_kde_rec_fftw_n-in)
figure(5)
plot(in_kde_rec_n,'-r')
hold on
plot(in,'-g')
plot(in_kde_rec_fftw_n,'-b')
hold off

figure(6)
plot(in_kde_rec,'-r')
hold on
plot(in_kde_rec_fftw,'-b')
hold off

N=length(initial_data);
y1=dct(initial_data*sqrt(N));
y2=kde_dct1d(initial_data);
figure(1)
plot(y1,'-r')
hold on
plot(y2,'-b')
hold off

%%
dct_error.y=abs(y1-y2)./abs(y1);
dct_error.mean=mean(dct_error.y);

dct_error.max=max(dct_error.y);

dct_error

figure(2)
plot(dct_error.y)
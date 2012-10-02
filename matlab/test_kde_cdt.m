N=length(initial_data);
y1=dct(initial_data);
y2=kde_dct1d(initial_data)/sqrt(N);
figure
plot(y1,'-r')
hold on
plot(y2,'-b')
hold off
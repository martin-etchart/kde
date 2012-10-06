clear all; close all; clc;

if 0
	randn1 = randn(100,1);
	randn2 = randn(100,1);
	randn3 = randn(100,1);
	
	data=[randn1;randn2*2+35 ;randn3+55];
	
	save('data','data')
	
	fid = fopen('data.txt','w'); fprintf(fid,'%2.16f\n',data);
	fclose(fid);
else
	load('data.mat');
end

[h,x]=hist(data,256);
dx=x(2)-x(1);
h=h/sum(h*dx);
[bw,f,xi,v_cdf]=kde2(data,128,min(data)-5,max(data)+5);


figure(1)
bar(x,h);
hold on
plot(xi,f,'-r');
hold off

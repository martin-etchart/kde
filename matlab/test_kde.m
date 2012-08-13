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

kde(data,2^14,min(data)-5,max(data)+5);

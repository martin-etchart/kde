close all
clear all
clc

% filename = 'trimodal2gaussian.png';
filename = 'peppers.png';

I=im2double(rgb2gray(imread(filename)));

% probando guardar y levantar imagenes como txt
jfc_vector_save_simple_mat(I,[filename(1:end-3) 'txt']);
Ir = jfc_vector_read_simple_mat([filename(1:end-3) 'txt']);

% probando el histograma de c
cBins=load('cBins.txt');
cHist=load('cHist.txt');
figure; stem(cBins,cHist); title('Histograma en C')

% [matlabHist,matlabBins] = imhist(I,11);
% dif = cHist'-matlabHist;
% fprintf('La diferencia del histograma da: %f\n',sum(dif));
% disp(dif')
% disp('      C           Matlab')
% disp([cHist' matlabHist])

[Iseg3,sep3] = otsu(I,3);
[Iseg2,sep2] = otsu(I,2);
figure; imshow([I Iseg3;Iseg2 zeros(size(I))]); 
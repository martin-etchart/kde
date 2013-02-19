close all
clear all
clc

filename = '../pics/trimodal2gaussian.png';
% filename = 'peppers.png';

I=im2double(rgb2gray(imread(filename)));

% probando guardar y levantar imagenes como txt
jfc_vector_save_simple_mat(I,[filename(1:end-3) 'txt']);
Ir = jfc_vector_read_simple_mat([filename(1:end-3) 'txt']);

% probando el histograma de c
cBins=load('../c/@build/cBins.txt');
cHist=load('../c/@build/cHist.txt');
figure; stem(cBins,cHist,'r'); title('Histograma en C')

% [matlabHist,matlabBins] = imhist(I,11);
% dif = cHist'-matlabHist;
% fprintf('La diferencia del histograma da: %f\n',sum(dif));
% disp(dif')
% disp('      C           Matlab')
% disp([cHist' matlabHist])

PC=load('../c/@build/P.txt');
wC=load('../c/@build/w.txt');
muC=load('../c/@build/mu.txt');
figure;plot(wC);title('Funcion de probabilidad')
sigma2BC=load('../c/@build/sigma2B.txt');
figure;plot(sigma2BC);title('sigma2B en C')

[Iseg3,sep3] = otsu(I,2);
% [Iseg2,sep2] = otsu(I,2);
% figure; imshow([I Iseg3;Iseg2 zeros(size(I))]);
figure; imshow([I Iseg3]); 
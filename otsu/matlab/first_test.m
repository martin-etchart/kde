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
nbins=length(cBins);
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

w0Caux=load('../c/@build/w0.txt');
w2Caux=load('../c/@build/w2.txt');
w0C=zeros(nbins,nbins);
w2C=zeros(nbins,nbins);
for i=1:size(w0C,1)
    w0C(i,:)=w0Caux((i-1)*size(w0C,2)+1:i*size(w0C,2));
    w2C(i,:)=w2Caux((i-1)*size(w2C,2)+1:i*size(w2C,2));
end

mu0Caux=load('../c/@build/mu0.txt');
mu2Caux=load('../c/@build/mu2.txt');
mu0C=zeros(nbins,nbins);
mu2C=zeros(nbins,nbins);
for i=1:size(mu0C,1)
    mu0C(i,:)=mu0Caux((i-1)*size(mu0C,2)+1:i*size(mu0C,2));
    mu2C(i,:)=mu2Caux((i-1)*size(mu2C,2)+1:i*size(mu2C,2));
end






data_outC=load('../c/@build/data_out.txt');
IsegC=zeros(size(I));
for i=1:size(I,1)
    IsegC(i,:)=data_outC((i-1)*size(I,2)+1:i*size(I,2));
end

[Iseg,sep] = otsu(I,3);
% figure; imshow([I Iseg;IsegC zeros(size(I))]);
figure('Name', 'Imagen Original - Imagen segmentada en Matlab - Imagen segmentada en C'); imshow([I Iseg IsegC]);
% figure; imshow([I Iseg3]); 

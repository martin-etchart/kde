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


data_outC=load('../c/@build/data_out.txt');
IsegC=zeros(size(I));
for i=1:size(I,1)
    IsegC(i,:)=data_outC((i-1)*size(I,2)+1:i*size(I,2));
end

[Iseg,sep] = otsu(I,2);
% figure; imshow([I Iseg;IsegC zeros(size(I))]);
figure('Name', 'Imagen Original - Imagen segmentada en Matlab - Imagen segmentada en C'); imshow([I Iseg IsegC]);
% figure; imshow([I Iseg3]); 

close all
clear all
clc

filename = 'trimodal.png';
% filename = 'peppers.png';

I=im2double(rgb2gray(imread(filename)));

% probando guardar y levantar imagenes como txt
jfc_vector_save_simple_mat(I,[filename(1:end-3) 'txt']);
Ir = jfc_vector_read_simple_mat([filename(1:end-3) 'txt']);
figure; imshow([I Ir]); 

% probando el histograma de c
% cBins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
% cHist = [183643, 0, 0, 0, 0, 22198, 0, 0, 0, 0, 56303];
cBins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
cHist = [200, 10718, 92024, 28652, 24083, 17641, 9565, 6324, 3456, 2445, 1500];
[matlabHist,matlabBins] = imhist(I,11);
dif = cHist'-matlabHist;
fprintf('La diferencia del histograma da: %f\n',sum(dif));
disp(dif')
disp('      C           Matlab')
disp([cHist' matlabHist])

close all
clear all
home

% the picture scales to fill the whole window
iptsetpref('ImshowBorder','tight');
% removes menu and toolbar from all new figures
% set(0,'DefaultFigureMenu','none');
%makes disp() calls show things without empty lines
format compact;

addpath(genpath('/home/mat/develop/matlab/oct-mat/tmaps/'))
addpath(genpath('/home/mat/develop/kde/'))

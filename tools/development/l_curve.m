
clear all
close all
clc


regu = [1e-2;  1e-4;  1e-6;  1e-8];
misfit = [1.355811e+01;  2.799345e+00;  4.047851e+00;  4.099760e+00];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


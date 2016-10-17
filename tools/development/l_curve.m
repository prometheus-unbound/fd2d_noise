
clear all
close all
clc

%% random 0.10
% structure cc
regu = [1e-1;  1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [1.093291e+02; 2.226342e+01; 6.193898e+00; 5.039486e+00; 5.329827e+00; 4.180563e+00];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


% source ampdiff equal cc
regu = [1e-1;  1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [1.900206e+01; 1.881410e+01; 1.755263e+01; 1.181557e+01; 3.567502e+00; 1.256639e+00];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


% source loga equal homogeneous
regu = [1e-1;  1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [6.387155e+01; 5.631987e+01; 3.219082e+01; 9.807957e+00; 2.502154e+00; 7.205383e-01];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


% source loga equal cc
regu = [1e-1;  1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [6.307025e+01; 5.515670e+01; 3.031675e+01; 8.325221e+00; 1.583784e+00; 3.238499e-01];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');





%% random 0.07
% random 0.07 - structure cc
regu = [1e-1;  1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [2.311953e+02; 7.149076e+01; 5.257688e+01; 4.885496e+01; 5.088960e+01; 5.118338e+01];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


% source ampdiff equal cc
regu = [1e-1;  1e-2;  1e-3;  1e-4;  1e-5;];
misfit = [2.555731e+02; 1.756853e+02; 8.972398e+01; 5.310749e+01; 3.206296e+01;];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


% source loga equal cc
regu = [1e-1;  1e-2;  1e-3;  1e-4;  1e-5;];
misfit = [5.938866e+01; 5.265943e+01; 3.147307e+01; 1.078507e+01; 3.551920e+00;];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


% source loga equal homogeneous
regu = [1e-1;  1e-2;  1e-3;  1e-4;  1e-5;];
misfit = [6.046780e+01; 5.391802e+01; 3.270924e+01; 1.222195e+01; 4.705168e+00;];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');




%% test
regu_source =    [1e-1;  1e-5;  1e-5;  1e-5;  1e-6;  1e-6;  1e-6;];
regu_structure = [1e-1;  1e-2;  1e-3;  1e-4;  1e-2;  1e-3;  1e-4;];
misfit = [1.243704e+02;  3.212904e+01;  1.718684e+01;  7.489421e+00;  3.173638e+01;  1.598900e+01;  6.491985e+00];

plot3( log10(regu_source), log10(regu_structure), misfit, 'x' )




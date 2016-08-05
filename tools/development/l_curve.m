
clear all
close all
clc

% smooth loga
regu = [1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [8.568545e+01;  4.586514e+01;  1.137602e+01;  1.969683e+00;  4.375895e-01];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


% smooth ampdiff
regu = [1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [2.936537e+01;  2.699236e+01;  1.793781e+01;  5.668404e+00;  1.908705e+00];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');





% coarse loga
regu = [1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [8.163603e+01;  3.873114e+01;  8.436142e+00;  1.356098e+00;  2.594611e-01];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');


% coarse ampdiff
regu = [1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [2.924174e+01;  2.631791e+01;  1.640061e+01;  4.352361e+00;  1.019265e+00];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');





% smooth wd
regu = [1e-2;  1e-3;  1e-4;  1e-5;  1e-6];
misfit = [1.310099e+02;  1.282889e+02;  1.157190e+02;  1.005077e+02;  8.944697e+01];

figure
semilogx(regu, misfit,'*')
set(gca,'XDir','Reverse');
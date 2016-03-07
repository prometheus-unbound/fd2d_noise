
clear all
close all
clc


% traditional
regu_t = [1e-0; 1e-1; 1e-2; 5e-1; 5e-2];
misfit_t = [5.867381e+01; 4.706844e+01; 4.087419e+01; 5.556599e+01; 4.434295e+01];

figure
semilogx(regu_t, misfit_t,'*')
set(gca,'XDir','Reverse');



% andreas idea
regu_a = [1e-0; 1e-1; 1e-2; 5e-1; 5e-2];
misfit_a = [2.633027e+01; 1.405917e+01; 8.185111e+00; 2.274586e+01; 1.146986e+01];

figure
semilogx(regu_a, misfit_a,'*')
set(gca,'XDir','Reverse');



% source point1
regu_a = [1e-4; 1e-5; 1e-6];
misfit_a = [2.681341e+00; 3.626467e-01; 5.369712e-02];

figure
semilogx(regu_a, misfit_a,'*')
set(gca,'XDir','Reverse');



regu_a = [1e-0; 1e-3; 1e-4];
misfit_a = [1.129963e+02; 3.420132e+01; 1.608990e+01];

figure
semilogx(regu_a, misfit_a,'*')
set(gca,'XDir','Reverse');


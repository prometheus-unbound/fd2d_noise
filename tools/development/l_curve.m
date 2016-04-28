
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
regu_a = [1e-0; 1e-2; 1e-4];
misfit_a = [9.727298e-01; 9.558338e-01; 6.656615e-01];

figure
semilogx(regu_a, misfit_a,'*')
set(gca,'XDir','Reverse');



regu_a = [1e-0; 1e-1; 1e-2; 7e-3; 5e-3; 3e-3; 1e-3];
misfit_a = [1.244949e+00; 6.059158e-01; 3.298011e-01; 3.127824e-01; 2.994536e-01; 2.824631e-01; 2.506788e-01;];

figure
semilogx(regu_a, misfit_a,'*')
set(gca,'XDir','Reverse');


regu_cc = [1e-0; 1e-1; 1e-2; 1e-3];
misfit_cc = [5.364138e+00; 2.122951e+00; 1.702147e+00; 1.540107e+00];

figure
semilogx(regu_cc, misfit_cc,'*')
set(gca,'XDir','Reverse');




regu_ampdiffsou = [1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7];
misfit_ampdiffsou = [2.112548e+01; 1.969718e+01; 1.339247e+01; 3.198018e+00; 4.408981e-01; 6.458783e-02];

figure
semilogx(regu_ampdiffsou, misfit_ampdiffsou,'*')
set(gca,'XDir','Reverse');



regu_ampdiffsou = [1e+2; 1e+1; 1e0; 1e-1; 1e-2; 1e-3; 1e-4];
misfit_ampdiffsou = [2.668754e+00; 2.587642e+00; 2.116070e+00; 1.246788e+00; 6.884613e-01; 4.171481e-01; 2.719252e-01];

figure
semilogx(regu_ampdiffsou, misfit_ampdiffsou,'*')
set(gca,'XDir','Reverse');



regu_ampdiffsou = [1e+2; 1e+1; 1e0; 1e-1; 1e-2; 1e-3];
misfit_ampdiffsou = [1.361341e+01; 1.271753e+01; 5.364138e+00; 2.122951e+00; 1.702147e+00; 1.540107e+00];

figure
semilogx(regu_ampdiffsou, misfit_ampdiffsou,'*')
set(gca,'XDir','Reverse');



regu_ampdiffsou = [1e0; 1e-2; 1e-4];
misfit_ampdiffsou = [9.727298e-01; 9.558338e-01; 6.656411e-01];

figure
semilogx(regu_ampdiffsou, misfit_ampdiffsou,'*')
set(gca,'XDir','Reverse');


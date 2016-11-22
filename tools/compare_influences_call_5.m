
clear all


folder = '~/Desktop/runs/paper/data_inversion/';
random = [0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.20];


% load reference
equal_homog = load([folder 'data_16_ref_0_equal_homogeneous.mat']);
gaussian_homog = load([folder 'data_16_ref_0_gaussian_homogeneous.mat']);
t = equal_homog.t;
load([folder 'array_16_ref.mat']);


% load data files
for i = 1:length(random)
    
    equal( i ) = load([folder 'data_16_ref_0_equal_random_' sprintf('%3.2f', random(i)) '_0.8e10_nosmooth.mat'], 'c_data');
    gaussian( i ) = load([folder 'data_16_ref_0_gaussian_random_' sprintf('%3.2f', random(i)) '_0.8e10_nosmooth.mat'], 'c_data');
    
end


% compute misfits
misfit_eqho_vs_eq = zeros(length(random), 4);
misfit_eqho_vs_gauss = zeros(length(random), 4);
misfit_gaussho_vs_gauss = zeros(length(random), 4);
misfit_eq_vs_gauss = zeros(length(random), 4);

for i = 1:length(random)
   
    misfit_eqho_vs_eq(i,:) = sum( compare_influences( equal_homog.c_data, equal(i).c_data, t, array, ref_stat, 'no' ), 1);
    misfit_eqho_vs_gauss(i,:) = sum( compare_influences( equal_homog.c_data, gaussian(i).c_data, t, array, ref_stat, 'no' ), 1);
    misfit_gaussho_vs_gauss(i,:) = sum( compare_influences( gaussian_homog.c_data, gaussian(i).c_data, t, array, ref_stat, 'no' ), 1);
    misfit_eq_vs_gauss(i,:) = sum( compare_influences( equal(i).c_data, gaussian(i).c_data, t, array, ref_stat, 'no' ), 1);
    
end


% plot misfits
titles(1,:) = 'loga';
titles(2,:) = 'amp ';
titles(3,:) = 'cc  ';
titles(4,:) = 'wd  ';

indices = 4:length(random);
for i = 1:4
    
    subplot( 2, 2, i )
    hold on
    plot( random(indices), misfit_eqho_vs_eq(indices,i), 'r' )
    plot( random(indices), misfit_eqho_vs_gauss(indices,i), 'b' )
%     plot( random(indices), misfit_gaussho_vs_gauss(indices,i), 'k' )
%     plot( random(indices), misfit_eq_vs_gauss(indices,i), 'c' )
    grid on
    
    title( titles(i,:), 'Interpreter', 'none' )
    
end

legend( 'homog. structure, equal noise dist.  VS  random structure, equal noise dist.', 'homog. structure, equal noise dist.  VS  random structure, gaussian noise dist.' );% , 'Location', 'Outside' )

return


path = '~/Desktop/model_25.mat';
load(path)
misfit_test1(1,:) = sum( compare_influences( equal_homog.c_data, gaussian_homog.c_data, t, array, ref_stat, 'no' ), 1);
misfit_test2(1,:) = sum( compare_influences( gaussian(4).c_data, model.correlation, t, array, ref_stat, 'no' ), 1);
misfit_test3(1,:) = sum( compare_influences( gaussian_homog.c_data, model.correlation, t, array, ref_stat, 'no' ), 1);

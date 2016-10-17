
clear all


folder = '~/Desktop/runs/paper/data/';
random = [0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.20];


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
misfit_eqhovseq = zeros(length(random), 4);
misfit_eqhovsgauss = zeros(length(random), 4);
misfit_gausshovsgauss = zeros(length(random), 4);
misfit_eqvsgauss = zeros(length(random), 4);

for i = 1:length(random)
   
    misfit_eqhovseq(i,:) = sum( compare_influences( equal_homog.c_data, equal(i).c_data, t, array, ref_stat, 'no' ), 1);
    misfit_eqhovsgauss(i,:) = sum( compare_influences( equal_homog.c_data, gaussian(i).c_data, t, array, ref_stat, 'no' ), 1);
    misfit_gausshovsgauss(i,:) = sum( compare_influences( gaussian_homog.c_data, gaussian(i).c_data, t, array, ref_stat, 'no' ), 1);
    misfit_eqvsgauss(i,:) = sum( compare_influences( equal(i).c_data, gaussian(i).c_data, t, array, ref_stat, 'no' ), 1);
    
end


% plot misfits
titles(1,:) = 'loga';
titles(2,:) = 'amp ';
titles(3,:) = 'cc  ';
titles(4,:) = 'wd  ';

for i = 1:4
    
    subplot( 2, 2, i )
    hold on
    plot( random, misfit_eqhovseq(:,i), 'r' )
    plot( random, misfit_eqhovsgauss(:,i), 'b' )
    plot( random, misfit_gausshovsgauss(:,i), 'k' )
    plot( random, misfit_eqvsgauss(:,i), 'c' )
    grid on
    
    title( titles(i,:), 'Interpreter', 'none' )
    
end




path = '~/Desktop/model_25.mat';
load(path)
misfit_test1(1,:) = sum( compare_influences( equal_homog.c_data, gaussian_homog.c_data, t, array, ref_stat, 'no' ), 1);
misfit_test2(1,:) = sum( compare_influences( gaussian(4).c_data, model.correlation, t, array, ref_stat, 'no' ), 1);
misfit_test3(1,:) = sum( compare_influences( gaussian_homog.c_data, model.correlation, t, array, ref_stat, 'no' ), 1);


clear all
clc


u_0_homog_h = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_homog.mat');
u_0_homog_h1g = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_homog.mat'); 
% u_0_random_h_0001 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.001_norm.mat');
% u_0_random_h1g_0001 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.001_norm.mat');
% u_0_random_h_001 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.01_norm.mat');
% u_0_random_h1g_001 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.01_norm.mat');
u_0_random_h_005 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.05_norm.mat');
u_0_random_h1g_005 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.05_norm.mat');
u_0_random_h_01 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.1_norm.mat');
u_0_random_h1g_01 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.1_norm.mat');
u_0_random_h_025 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.25_norm.mat');
u_0_random_h1g_025 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.25_norm.mat');
u_0_random_h_02 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.2_norm.mat');
u_0_random_h1g_02 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.2_norm.mat');
u_0_random_h_03 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.3_norm.mat');
u_0_random_h1g_03 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.3_norm.mat');
u_0_random_h_04 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.4_norm.mat');
u_0_random_h1g_04 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.4_norm.mat');
u_0_random_h_05 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h_random_0.5_norm.mat');
u_0_random_h1g_05 = load('~/Desktop/runs/2016_start/data_random_3/data_16_ref_0_1h1g_random_0.5_norm.mat');

t = u_0_random_h_005.t;
load('~/Desktop/runs/2016_start/data_iugg/array_16_ref.mat');

% misfit_0001 = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_0001.c_data, t, array, ref_stat, 'no' );
% misfit_001 = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_001.c_data, t, array, ref_stat, 'no' );
misfit_005 = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_005.c_data, t, array, ref_stat, 'yes' );
misfit_01  = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_01.c_data, t, array, ref_stat, 'no' );
misfit_02  = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_02.c_data, t, array, ref_stat, 'no' );
misfit_025  = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_025.c_data, t, array, ref_stat, 'no' );
misfit_03  = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_03.c_data, t, array, ref_stat, 'no' );
misfit_04  = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_04.c_data, t, array, ref_stat, 'no' );
misfit_05  = compare_influences( u_0_homog_h.c_data, u_0_random_h1g_05.c_data, t, array, ref_stat, 'no' );


titles(1,:) = 'loga';
titles(2,:) = 'amp ';
titles(3,:) = 'cc  ';
titles(4,:) = 'wd  ';

color = 'b*-';
figure(2)
for i=1:4
    subplot(2,2,i)
    hold on
    handle(1,:) = plot([0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5], [sum(misfit_005(:,i)), sum(misfit_01(:,i)), sum(misfit_02(:,i)), sum(misfit_025(:,i)), sum(misfit_03(:,i)), sum(misfit_04(:,i)), sum(misfit_05(:,i))], color);
    title(titles(i,:))
end




% misfit_0001_2 = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_0001.c_data, t, array, ref_stat, 'no' );
% misfit_001_2 = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_001.c_data, t, array, ref_stat, 'no' );
misfit_005_2 = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_005.c_data, t, array, ref_stat, 'no' );
misfit_01_2  = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_01.c_data, t, array, ref_stat, 'no' );
misfit_02_2  = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_02.c_data, t, array, ref_stat, 'no' );
misfit_025_2  = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_025.c_data, t, array, ref_stat, 'no' );
misfit_03_2  = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_03.c_data, t, array, ref_stat, 'no' );
misfit_04_2  = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_04.c_data, t, array, ref_stat, 'no' );
misfit_05_2  = compare_influences( u_0_homog_h1g.c_data, u_0_random_h1g_05.c_data, t, array, ref_stat, 'no' );


color = 'r+-';
figure(2)
for i=1:4
    subplot(2,2,i)
    hold on
    handle(2,:) = plot([0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5], [sum(misfit_005_2(:,i)), sum(misfit_01_2(:,i)), sum(misfit_02_2(:,i)), sum(misfit_025_2(:,i)), sum(misfit_03_2(:,i)), sum(misfit_04_2(:,i)), sum(misfit_05_2(:,i))], color);
    title(titles(i,:))
    
end

legend(handle, 'homog h to structure h1g', 'homog h1g to structure h1g')
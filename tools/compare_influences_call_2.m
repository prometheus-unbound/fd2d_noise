
clear all
clc


u_0_homog_h = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_0.3_1_homog.mat');
u_0_homog_h2g = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_0.3_1_homog.mat'); 
u_0_random_h_0001 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.001.mat');
u_0_random_h2g_0001 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.001.mat');
u_0_random_h_001 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.01.mat');
u_0_random_h2g_001 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.01.mat');
u_0_random_h_005 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.05.mat');
u_0_random_h2g_005 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.05.mat');
u_0_random_h_01 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.1.mat');
u_0_random_h2g_01 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.1.mat');
u_0_random_h_025 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.25.mat');
u_0_random_h2g_025 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.25.mat');
u_0_random_h_02 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.2.mat');
u_0_random_h2g_02 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.2.mat');
u_0_random_h_03 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.3.mat');
u_0_random_h2g_03 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.3.mat');
u_0_random_h_04 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.4.mat');
u_0_random_h2g_04 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.4.mat');
u_0_random_h_05 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h_random_0.5.mat');
u_0_random_h2g_05 = load('~/Desktop/runs/2016_start/data_random_h2g/data_16_ref_91_h2g_random_0.5.mat');

t = u_0_random_h_005.t;
load('~/Desktop/runs/2016_start/data_iugg/array_16_ref.mat');

misfit_0001 = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_0001.c_data, t, array, ref_stat, 'no' );
misfit_001 = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_001.c_data, t, array, ref_stat, 'no' );
misfit_005 = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_005.c_data, t, array, ref_stat, 'yes' );
misfit_01  = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_01.c_data, t, array, ref_stat, 'no' );
misfit_02  = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_02.c_data, t, array, ref_stat, 'no' );
misfit_025  = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_025.c_data, t, array, ref_stat, 'no' );
misfit_03  = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_03.c_data, t, array, ref_stat, 'no' );
misfit_04  = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_04.c_data, t, array, ref_stat, 'no' );
misfit_05  = compare_influences( u_0_homog_h.c_data, u_0_random_h2g_05.c_data, t, array, ref_stat, 'no' );


titles(1,:) = 'loga';
titles(2,:) = 'amp ';
titles(3,:) = 'cc  ';
titles(4,:) = 'wd  ';

color = 'b*-';
figure(1)
for i=1:4
    subplot(2,2,i)
    hold on
    handle(1,:) = plot([0.001, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5], [sum(misfit_0001(:,i)), sum(misfit_001(:,i)), sum(misfit_005(:,i)), sum(misfit_01(:,i)), sum(misfit_02(:,i)), sum(misfit_025(:,i)), sum(misfit_03(:,i)), sum(misfit_04(:,i)), sum(misfit_05(:,i))], color);
    title(titles(i,:))
end




misfit_0001_2 = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_0001.c_data, t, array, ref_stat, 'no' );
misfit_001_2 = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_001.c_data, t, array, ref_stat, 'no' );
misfit_005_2 = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_005.c_data, t, array, ref_stat, 'yes' );
misfit_01_2  = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_01.c_data, t, array, ref_stat, 'no' );
misfit_02_2  = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_02.c_data, t, array, ref_stat, 'no' );
misfit_025_2  = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_025.c_data, t, array, ref_stat, 'no' );
misfit_03_2  = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_03.c_data, t, array, ref_stat, 'no' );
misfit_04_2  = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_04.c_data, t, array, ref_stat, 'no' );
misfit_05_2  = compare_influences( u_0_homog_h2g.c_data, u_0_random_h2g_05.c_data, t, array, ref_stat, 'no' );


color = 'r+-';
figure(1)
for i=1:4
    subplot(2,2,i)
    hold on
    handle(2,:) = plot([0.001, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5], [sum(misfit_0001_2(:,i)), sum(misfit_001_2(:,i)), sum(misfit_005_2(:,i)), sum(misfit_01_2(:,i)), sum(misfit_02_2(:,i)), sum(misfit_025_2(:,i)), sum(misfit_03_2(:,i)), sum(misfit_04_2(:,i)), sum(misfit_05_2(:,i))], color);
    title(titles(i,:))
    
end

legend(handle, 'homog h to structure h2g', 'homog h2g to structure h2g')
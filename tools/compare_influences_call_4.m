
clear all
% clc


u_homog_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_homog.mat');
u_homog_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_homog.mat');



u_10_iaf_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_10_in_array_few.mat'); 
u_10_iam_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_10_in_array_many.mat'); 
u_10_ltaf_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_10_lt_of_array_few.mat'); 
u_10_ltam_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_10_lt_of_array_many.mat');

u_20_iaf_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_20_in_array_few.mat'); 
u_20_iam_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_20_in_array_many.mat'); 

u_15_iaf_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_15_in_array_few.mat'); 
u_15_iam_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_15_in_array_many.mat'); 

u_5_iaf_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_5_in_array_few.mat'); 
u_5_iam_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_5_in_array_many.mat'); 

u_3_iaf_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_3_in_array_few.mat'); 
u_3_iam_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_3_in_array_many.mat'); 



u_10_iaf_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_10_in_array_few.mat'); 
u_10_iam_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_10_in_array_many.mat'); 
u_10_ltaf_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_10_lt_of_array_few.mat'); 
u_10_ltam_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_10_lt_of_array_many.mat'); 

u_20_iaf_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_20_in_array_few.mat'); 
u_20_iam_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_20_in_array_many.mat'); 

u_15_iaf_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_15_in_array_few.mat'); 
u_15_iam_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_15_in_array_many.mat'); 

u_5_iaf_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_5_in_array_few.mat'); 
u_5_iam_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_5_in_array_many.mat'); 

u_3_iaf_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_3_in_array_few.mat'); 
u_3_iam_h1g = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h1g_random_3_in_array_many.mat'); 



u_10_10_iaf_h = load('~/Desktop/runs/2016_start/data_scatterer_10e9/data_16_ref_0_h_random_10_in_array_few.mat'); 
u_10_10_iam_h = load('~/Desktop/runs/2016_start/data_scatterer_10e9/data_16_ref_0_h_random_10_in_array_many.mat'); 
u_10_10_ltaf_h = load('~/Desktop/runs/2016_start/data_scatterer_10e9/data_16_ref_0_h_random_10_lt_of_array_few.mat'); 
u_10_10_ltam_h = load('~/Desktop/runs/2016_start/data_scatterer_10e9/data_16_ref_0_h_random_10_lt_of_array_many.mat');

u_10_10_iaf_h1g = load('~/Desktop/runs/2016_start/data_scatterer_10e9/data_16_ref_0_h1g_random_10_in_array_few.mat'); 
u_10_10_iam_h1g = load('~/Desktop/runs/2016_start/data_scatterer_10e9/data_16_ref_0_h1g_random_10_in_array_many.mat'); 
u_10_10_ltaf_h1g = load('~/Desktop/runs/2016_start/data_scatterer_10e9/data_16_ref_0_h1g_random_10_lt_of_array_few.mat'); 
u_10_10_ltam_h1g = load('~/Desktop/runs/2016_start/data_scatterer_10e9/data_16_ref_0_h1g_random_10_lt_of_array_many.mat'); 

t = u_homog_h.t;
load('~/Desktop/runs/2016_start/data_iugg/array_16_ref.mat');


misfit_3_iam = compare_influences( u_homog_h1g.c_data, u_3_iam_h1g.c_data, t, array, ref_stat, 'no' );
misfit_5_iam = compare_influences( u_homog_h1g.c_data, u_5_iam_h1g.c_data, t, array, ref_stat, 'no' );
misfit_10_iam = compare_influences( u_homog_h1g.c_data, u_10_iam_h1g.c_data, t, array, ref_stat, 'no' );
misfit_15_iam = compare_influences( u_homog_h1g.c_data, u_15_iam_h1g.c_data, t, array, ref_stat, 'no' );
misfit_20_iam = compare_influences( u_homog_h1g.c_data, u_20_iam_h1g.c_data, t, array, ref_stat, 'no' );

misfit_3_iaf = compare_influences( u_homog_h1g.c_data, u_3_iaf_h1g.c_data, t, array, ref_stat, 'no' );
misfit_5_iaf = compare_influences( u_homog_h1g.c_data, u_5_iaf_h1g.c_data, t, array, ref_stat, 'no' );
misfit_10_iaf = compare_influences( u_homog_h1g.c_data, u_10_iaf_h1g.c_data, t, array, ref_stat, 'no' );
misfit_15_iaf = compare_influences( u_homog_h1g.c_data, u_15_iaf_h1g.c_data, t, array, ref_stat, 'no' );
misfit_20_iaf = compare_influences( u_homog_h1g.c_data, u_20_iaf_h1g.c_data, t, array, ref_stat, 'no' );

misfit_10_ltam = compare_influences( u_homog_h1g.c_data, u_10_ltam_h1g.c_data, t, array, ref_stat, 'no' );
misfit_10_ltaf = compare_influences( u_homog_h1g.c_data, u_10_ltaf_h1g.c_data, t, array, ref_stat, 'no' );


misfit_10_10 = compare_influences( u_homog_h1g.c_data, u_10_10_iam_h1g.c_data, t, array, ref_stat, 'no' );


titles(1,:) = 'loga';
titles(2,:) = 'amp ';
titles(3,:) = 'cc  ';
titles(4,:) = 'wd  ';

color = 'b*-';
figure(1)
for i=1:4
    subplot(2,2,i)
    hold on
    handle(1,:) = plot([3 5 10 15 20], [sum(misfit_3_iam(:,i)), sum(misfit_5_iam(:,i)), sum(misfit_10_iam(:,i)), sum(misfit_15_iam(:,i)), sum(misfit_20_iam(:,i))], 'r*-');
    handle(2,:) = plot([3 5 10 15 20], [sum(misfit_3_iaf(:,i)), sum(misfit_5_iaf(:,i)), sum(misfit_10_iaf(:,i)), sum(misfit_15_iaf(:,i)), sum(misfit_20_iaf(:,i))], 'b*-');
    
    handle(4,:) = plot([10], [sum(misfit_10_ltaf(:,i))], 'y*-');
    handle(3,:) = plot([10], [sum(misfit_10_ltam(:,i))], 'g*-');
    title(titles(i,:))
end






misfit_3_iam = compare_influences( u_homog_h.c_data, u_3_iam_h.c_data, t, array, ref_stat, 'no' );
misfit_5_iam = compare_influences( u_homog_h.c_data, u_5_iam_h.c_data, t, array, ref_stat, 'no' );
misfit_10_iam = compare_influences( u_homog_h.c_data, u_10_iam_h.c_data, t, array, ref_stat, 'no' );
misfit_15_iam = compare_influences( u_homog_h.c_data, u_15_iam_h.c_data, t, array, ref_stat, 'no' );
misfit_20_iam = compare_influences( u_homog_h.c_data, u_20_iam_h.c_data, t, array, ref_stat, 'no' );

misfit_3_iaf = compare_influences( u_homog_h.c_data, u_3_iaf_h.c_data, t, array, ref_stat, 'no' );
misfit_5_iaf = compare_influences( u_homog_h.c_data, u_5_iaf_h.c_data, t, array, ref_stat, 'no' );
misfit_10_iaf = compare_influences( u_homog_h.c_data, u_10_iaf_h.c_data, t, array, ref_stat, 'no' );
misfit_15_iaf = compare_influences( u_homog_h.c_data, u_15_iaf_h.c_data, t, array, ref_stat, 'no' );
misfit_20_iaf = compare_influences( u_homog_h.c_data, u_20_iaf_h.c_data, t, array, ref_stat, 'no' );

misfit_10_ltam = compare_influences( u_homog_h.c_data, u_10_ltam_h.c_data, t, array, ref_stat, 'no' );
misfit_10_ltaf = compare_influences( u_homog_h.c_data, u_10_ltaf_h.c_data, t, array, ref_stat, 'no' );


misfit_10_10 = compare_influences( u_homog_h.c_data, u_10_10_iam_h.c_data, t, array, ref_stat, 'no' );


titles(1,:) = 'loga';
titles(2,:) = 'amp ';
titles(3,:) = 'cc  ';
titles(4,:) = 'wd  ';

color = 'b*-';
figure(3)
for i=1:4
    subplot(2,2,i)
    hold on
    handle(1,:) = plot([3 5 10 15 20], [sum(misfit_3_iam(:,i)), sum(misfit_5_iam(:,i)), sum(misfit_10_iam(:,i)), sum(misfit_15_iam(:,i)), sum(misfit_20_iam(:,i))], 'r*-');
    handle(2,:) = plot([3 5 10 15 20], [sum(misfit_3_iaf(:,i)), sum(misfit_5_iaf(:,i)), sum(misfit_10_iaf(:,i)), sum(misfit_15_iaf(:,i)), sum(misfit_20_iaf(:,i))], 'b*-');
    
    handle(4,:) = plot([10], [sum(misfit_10_ltaf(:,i))], 'y*-');
    handle(3,:) = plot([10], [sum(misfit_10_ltam(:,i))], 'g*-');
    title(titles(i,:))
end


% u_homog_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_homog.mat');
% u_5_ilta_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_5_in_lt_array.mat'); 
% u_7_ilta_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_7_in_lt_array.mat'); 
% u_10_ilta_h = load('~/Desktop/runs/2016_start/data_scatterer_5e9/data_16_ref_0_h_random_10_in_lt_array.mat'); 
% 
% t = u_homog_h.t;
% load('~/Desktop/runs/2016_start/data_iugg/array_16_ref.mat');
% 
% misfit_5_ilta = compare_influences( u_homog_h.c_data, u_5_ilta_h.c_data, t, array, ref_stat, 'no' );
% misfit_7_ilta = compare_influences( u_homog_h.c_data, u_7_ilta_h.c_data, t, array, ref_stat, 'no' );
% misfit_10_ilta = compare_influences( u_homog_h.c_data, u_10_ilta_h.c_data, t, array, ref_stat, 'no' );


clear all
% clc


h_homog_c = load('~/Desktop/runs/paper/source_inversions/20/data_16_ref_0_h_homog_coarse.mat');
h_homog_s = load('~/Desktop/runs/paper/source_inversions/20/data_16_ref_0_h_homog_smooth.mat');
h_iugg_c = load('~/Desktop/runs/paper/source_inversions/20/data_16_ref_0_h_iugg_coarse.mat');
h_iugg_s = load('~/Desktop/runs/paper/source_inversions/20/data_16_ref_0_h_iugg_smooth.mat');

h1g_homog_c = load('~/Desktop/runs/paper/source_inversions/20/data_16_ref_0_h1g_homog_coarse.mat');
h1g_homog_s = load('~/Desktop/runs/paper/source_inversions/20/data_16_ref_0_h1g_homog_smooth.mat');
h1g_iugg_c = load('~/Desktop/runs/paper/source_inversions/20/data_16_ref_0_h1g_iugg_coarse.mat');
h1g_iugg_s = load('~/Desktop/runs/paper/source_inversions/20/data_16_ref_0_h1g_iugg_smooth.mat');

h1g_iugg_c_10 = load('~/Desktop/runs/paper/source_inversions/10/data_16_ref_0_h1g_iugg_coarse.mat');
h1g_iugg_s_10 = load('~/Desktop/runs/paper/source_inversions/10/data_16_ref_0_h1g_iugg_smooth.mat');

h1g_iugg_c_15 = load('~/Desktop/runs/paper/source_inversions/15/data_16_ref_0_h1g_iugg_coarse.mat');
h1g_iugg_s_15 = load('~/Desktop/runs/paper/source_inversions/15/data_16_ref_0_h1g_iugg_smooth.mat');



t = h_homog_c.t;
load('~/Desktop/runs/paper/source_inversions/array_16_ref.mat');


misfit_1 = compare_influences( h1g_homog_c.c_data, h1g_homog_s.c_data, t, array, ref_stat, 'no' );

misfit_2 = compare_influences( h_homog_c.c_data, h1g_homog_c.c_data, t, array, ref_stat, 'no' );
misfit_3 = compare_influences( h_homog_c.c_data, h1g_iugg_c.c_data, t, array, ref_stat, 'no' );

misfit_4 = compare_influences( h_homog_s.c_data, h1g_homog_s.c_data, t, array, ref_stat, 'no' );
misfit_5 = compare_influences( h_homog_s.c_data, h1g_iugg_s.c_data, t, array, ref_stat, 'no' );


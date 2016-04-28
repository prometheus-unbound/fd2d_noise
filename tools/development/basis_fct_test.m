
clear all
close all
clc

basis(1) = load('~/Desktop/runs/basis_fct/data_16_ref_5_h2g_0.3_1_homog_nfft_5.mat');
basis(2) = load('~/Desktop/runs/basis_fct/data_16_ref_10_h2g_0.3_1_homog_nfft_5.mat');
basis(3) = load('~/Desktop/runs/basis_fct/data_16_ref_20_h2g_0.3_1_homog_nfft_5.mat');
basis(4) = load('~/Desktop/runs/basis_fct/data_16_ref_30_h2g_0.3_1_homog_nfft_5.mat');
basis(5) = load('~/Desktop/runs/basis_fct/data_16_ref_40_h2g_0.3_1_homog_nfft_5.mat');
basis(6) = load('~/Desktop/runs/basis_fct/data_16_ref_50_h2g_0.3_1_homog_nfft_5.mat');
basis(7) = load('~/Desktop/runs/basis_fct/data_16_ref_60_h2g_0.3_1_homog_nfft_5.mat');
basis(8) = load('~/Desktop/runs/basis_fct/data_16_ref_70_h2g_0.3_1_homog_nfft_5.mat');
basis(9) = load('~/Desktop/runs/basis_fct/data_16_ref_80_h2g_0.3_1_homog_nfft_5.mat');
basis(10) = load('~/Desktop/runs/basis_fct/data_16_ref_91_h2g_0.3_1_homog_nfft_5.mat');



to_plot = [5 6 7 8 9 10];

for i = 1:length(to_plot)
    
    h(i,:) = plot_recordings( basis(to_plot(i)).c_data(1,:), basis(to_plot(i)).t, 'vel', rand(1,3), true);  
    strings(i,:) = sprintf('%2s', num2str( to_plot(i)) );
    
end


legend(h,strings)
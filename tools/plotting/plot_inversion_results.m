
clear all
close all
clc

% get configuration
[~,~,nx,nz,dt,nt,~,~,~,n_basis_fct] = input_parameters();
if(n_basis_fct == 0)
    n_basis_fct = 1;
end


% path = '~/Desktop/runs/inversion_newest/source/2h2g_lr_nover_f1_bam/log_a/';
% path = '~/Desktop/runs/inversion_newest/source/2h2g_f1_bam_richtig/log_a/';
% path = '~/Desktop/runs/inversion_newest/coverage/49_full/';
% path = '~/Desktop/runs/inversion_newest/tradeoff/rand20/';
% path = '~/Desktop/runs/inversion_newest/point1/1h/wd_freqsamp5_nosmoo/';
path = '~/Desktop/';

% n_models = length( dir([path 'model_*']) );
n_models = 42;

model_final = load([path 'model_' num2str(n_models-1) '.mat']);
dist_inverted = reshape(model_final.model.m,nx,nz,n_basis_fct);
% dist_inverted = reshape(model_final.m,nx,nz,n_basis_fct);
% dist_inverted = 4000 * (1+dist_inverted);

% cm = cbrewer('div','RdBu',100,'PCHIP');
% m = max(max(max(dist_inverted)));
% clim = [-1*m 1*m];

% mine = fspecial('gaussian',[75 75], 30);
% mine = fspecial('gaussian',[50 50], 20);
% mine = fspecial('gaussian',[25 25], 10);
% dist_inverted = imfilter( dist_inverted, mine, 'replicate' );

[dist_true,~,clim] = make_noise_source('yes');
% load ../output/interferometry/array_49_ref_coverage.mat

if( ~exist('cm','var') )
    cm = [];
end

if( exist('clim','var') )
    plot_noise_sources(dist_inverted,[],cm,[clim(1) clim(2)]);
else
    plot_noise_sources(dist_inverted,[],cm,[]);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% plot correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform_1gaussian_homogeneous.mat');
% initial = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform1_homogeneous.mat');
% t = -(nt-1)*dt:dt:(nt-1)*dt;
% 
% figure
% plot_recordings(data.c_data,t,'vel','k',true);
% plot_recordings(initial.c_data,t,'vel','r',true);
% plot_recordings(model_final.cn,t,'vel','b',true);



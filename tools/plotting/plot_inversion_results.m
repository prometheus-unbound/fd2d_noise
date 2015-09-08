
clear all
close all
clc

% get configuration
[~,~,nx,nz,dt,nt,~,~,~,n_basis_fct] = input_parameters();
if(n_basis_fct == 0)
    n_basis_fct = 1;
end


% path = '~/Desktop/runs/inversion_basis_fct/source/d0_u2g_h_filtered/i10_u_h/loga_c/f2_lbfgs_unfiltered/';
% path = '~/Desktop/runs/inversion_basis_fct/source/d0_u2g_h_filtered/i10_u_h/loga_c/f2_lbfgs/';
% path = '~/Desktop/runs/inversion_basis_fct/source/d0_u2g_h/i10_u_h/loga_c/';

% n_models = length( dir([path 'model_*']) );
% model_final = load([path 'model_' num2str(n_models-1) '.mat']);
% dist_inverted = reshape(model_final.xn,nx,nz,n_basis_fct);


model_final = load('~/Desktop/model_1.mat');
dist_inverted = reshape(model_final.xn,nx,nz,n_basis_fct);

% cm = cbrewer('div','RdBu',100,'PCHIP');
% m = max(max(max(dist_inverted)));
% clim = [-m m];


[dist_true,~,clim] = make_noise_source('yes');

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



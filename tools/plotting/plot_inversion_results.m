
clear all
% close all
clc

% get configuration
[~,~,nx,nz,dt,nt,~,~,~,n_basis_fct] = input_parameters();
if(n_basis_fct == 0)
    n_basis_fct = 1;
end

usr_par.type = 'source';
usr_par.kernel.imfilter = fspecial('gaussian',[50 50], 20);

path = '~/Desktop/model_1.mat';


% path = sprintf('~/Desktop/models/model_%i.mat',i);

% path = '~/Desktop/runs/inversion/models/true_mu_structure_1.mat';
% path = '~/Desktop/runs/inversion/models/true_source_uniform_blob3.mat';

% path = '~/Desktop/runs/inversion/structure_1/uniform_blob3/structure/true_source/wd/';
% path = '~/Desktop/runs/inversion/structure_1/uniform_blob3/structure/homog_source/wd/';
% path = '~/Desktop/runs/inversion/structure_1/uniform_blob3/structure/source_from_log_a/wd/';

% path = '~/Desktop/runs/inversion/structure_1/uniform_blob3/source/homog_structure/unconstrained/log_a/';
% path = '~/Desktop/runs/inversion/structure_1/uniform_blob3/source/true_structure/log_a/';

% path = '~/Desktop/runs/inversion_newest/source/2h2g_lr_nover_f1_bam/log_a/';
% path = '~/Desktop/runs/inversion_newest/source/2h2g_f1_bam_richtig/log_a/';

% path = '~/Desktop/runs/inversion_newest/coverage/16_orig/';
% path = '~/Desktop/runs/inversion_newest/coverage/16_full/';
% path = '~/Desktop/runs/inversion_newest/coverage/49_full/';

% path = '~/Desktop/runs/inversion_newest/tradeoff/rand20/';
% path = '~/Desktop/runs/inversion_newest/point1/1h/wd_freqsamp5_nosmoo/';
% path = '~/Desktop/';


% n_models = length( dir([path 'model_*']) );


% model_final = load([path 'model_' num2str(n_models-1) '.mat']);
load(path);
m_parameters = reshape( map_m_to_parameters(model.m, usr_par), nx, nz, n_basis_fct );


% [dist_true,~,clim] = make_noise_source('yes');
% load ~/Desktop/runs/inversion_newest/data/coverage/array_16_ref_coverage_orig.mat
% load ~/Desktop/runs/inversion_newest/data/coverage/array_16_ref_coverage.mat
% load ~/Desktop/runs/inversion_newest/data/coverage/array_49_ref_coverage.mat



load ~/Desktop/runs/inversion_new_world/data/array_16_ref.mat



% cm = cbrewer('div','RdBu',100,'PCHIP');
% m = max(max(max(dist_inverted)));
% clim = [-1*m 1*m];
if( ~exist('cm','var') )
    cm = [];
end

if( exist('clim','var') )
    plot_models( m_parameters, array, cm, [clim(1) clim(2)] );
else
    plot_models( m_parameters, array, cm, [] );
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



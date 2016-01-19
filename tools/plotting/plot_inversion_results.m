
clear all
% close all
clc

% get configuration
[~,~,nx,nz,dt,nt,~,~,~,n_basis_fct] = input_parameters();
if(n_basis_fct == 0)
    n_basis_fct = 1;
end

% for i = 1:31

% path = '~/Desktop/mu_point4.mat';

path = '~/Desktop/model_4.mat';

% path = '~/Desktop/model_55_point2_equal_wd.mat';
% path = '~/Desktop/model_61_point2_logas_wd.mat';
% path = '~/Desktop/model_60_point2_trues_wd.mat';

% path = '~/Desktop/model_20_point2_equal_wd.mat';
% path = '~/Desktop/model_20_point2_logas_wd.mat';

% path = '~/Desktop/point_tests/point3/model_36_point3_trues_wd.mat';

% path = '~/Desktop/point_tests/point4/model_2_point4_trues_cc.mat';



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
model_final = load(path);
% model_initial = load('~/Desktop/model_0.mat');

dist_inverted = reshape(model_final.model.m,nx,nz,n_basis_fct);
% dist_inverted = reshape(model_final.model.m-model_initial.model.m,nx,nz,n_basis_fct);
% dist_inverted = reshape(model_final.mu,nx,nz,n_basis_fct);
% dist_inverted = reshape(model_final.m,nx,nz,n_basis_fct);
% dist_inverted = reshape(model_final.xn,nx,nz,n_basis_fct);
% dist_inverted = reshape(model_final.mu,nx,nz,n_basis_fct);


% clip_value = prctile( abs(reshape(dist_inverted,[],1)), 99.7 );
% [i,j] = find( abs(dist_inverted) >= clip_value );
% dist_inverted(i,j) = 0;

dist_inverted = 4.8e10 * (1+dist_inverted);


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
    plot_models(dist_inverted,array,cm,[clim(1) clim(2)]);
else
    plot_models(dist_inverted,array,cm,[]);
    plot_models(dist_inverted,array,cm,[4.75 4.85]*1e10);
end

pause(0.1)

% end


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



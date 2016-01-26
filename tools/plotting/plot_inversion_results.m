
clear all
% close all
clc

% get configuration
[~,~,nx,nz,dt,nt,~,model_type,~,n_basis_fct] = input_parameters();
if(n_basis_fct == 0)
    n_basis_fct = 1;
end


path = '~/Desktop/model_3.mat';
% path2 = '~/Desktop/source_iugg.mat';


% n_models = length( dir([path 'model_*']) );
% model_final = load([path 'model_' num2str(n_models-1) '.mat']);
load(path);

if( exist('path2', 'var') )
    model2 = load(path2);
    model.m = model.m - model2.model.m;
end

if( strcmp(path(end-11:end),'solution.mat') )
    model.m = mfinal;
else
    usr_par.kernel.imfilter = model.imfilter;
end

m_parameters = reshape( map_m_to_parameters(model.m, usr_par), nx, nz, n_basis_fct + 1);


% mu = define_material_parameters(nx,nz,model_type,'no');
% m_parameters(:,:,2) = m_parameters(:,:,2) - mu;


% load ~/Desktop/runs/inversion_new_world/data/array_16_ref.mat
if( ~exist('array', 'var') )
    array = [];
end

% cm = cbrewer('div','RdBu',100,'PCHIP');
% m = max(max(max(abs(m_parameters))));
% m = mean(mean(mean(abs(m_parameters))));
% clim = [0.999*m 1.001*m];
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



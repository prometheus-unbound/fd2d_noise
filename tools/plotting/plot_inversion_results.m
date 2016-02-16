
clear all
% close all
clc

[Lx, Lz, nx, nz, dt, nt, order, model_type, source_type, n_basis_fct, fw_nth] = input_parameters();

path = '~/Desktop/model_56.mat';
% path2 = '~/Desktop/model_69.mat';


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
    [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);
end


% get true source and material
if( usr_par.config.n_basis_fct == 0 )
    m_true = zeros( usr_par.config.nx, usr_par.config.nz, 2 );
else
    m_true = zeros( usr_par.config.nx, usr_par.config.nz, usr_par.config.n_basis_fct+1 );
end

% m_true(:,:,1:end-1) = make_noise_source( source_type, usr_par.config.n_basis_fct );
% m_true(:,:,end) = define_material_parameters( usr_par.config.nx, usr_par.config.nz, model_type );
% m_true = map_m_to_parameters( map_parameters_to_m( m_true, usr_par ) , usr_par );
% m_parameters = map_m_to_parameters(model.m, usr_par) - m_true;

m_parameters = map_m_to_parameters(model.m, usr_par);

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
% t = data.t;
% 
% figure
% plot_recordings(data.c_data,t,'vel','k',true);
% plot_recordings(initial.c_data,t,'vel','r',true);
% plot_recordings(model_final.cn,t,'vel','b',true);




clear all
% close all
clc

[Lx, Lz, nx, nz, dt, nt, order, model_type, source_type, n_basis_fct, fw_nth] = input_parameters();


% path = '~/Desktop/runs/2016_start/source_inversions/homog_0.mat';
% path2 = '~/Desktop/runs/2016_start/source_inversions/true_structure_0.mat';

path = '~/Desktop/model_0.mat';
% path2 = '~/Desktop/model_32.mat';

difftrue = 'no';
gradient = 'yes';



%% load model
load(path);



%% set usr_par; solution.mat contains usr_par
if( strcmp(path(end-11:end),'solution.mat') )
    model.m = mfinal;
else
    usr_par.network = []; usr_par.data = [];
    
    usr_par.ring.switch = 'no';
    usr_par.ring.x_center_ring = 1.0e6;
    usr_par.ring.z_center_ring = 1.0e6;
    usr_par.ring.radius = 6.4e5;
    usr_par.ring.thickness = 2.0e5;
    usr_par.ring.taper_strength = 70e8;
        
%     usr_par.ring.switch = 'no';
%     usr_par.ring.x_center_ring = 1.0e6;
%     usr_par.ring.z_center_ring = 1.0e6;
%     usr_par.ring.radius = 6.9e5;
%     usr_par.ring.thickness = 1.0e5;
%     usr_par.ring.taper_strength = 10e8;

    % usr_par.ring = model.ring;
   
    if( isfield(model.imfilter, 'source') )
        usr_par.kernel.imfilter.source = model.imfilter.source;
        usr_par.kernel.imfilter.structure = model.imfilter.structure;
    else
        usr_par.kernel.imfilter.source = model.imfilter;
        usr_par.kernel.imfilter.structure = model.imfilter;
        % usr_par.config.n_basis_fct = model.config.n_basis_fct;
        usr_par.config.n_basis_fct = 5;
    end
    
    
    [usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);
end


%% plot difference to true model
if( strcmp(difftrue,'yes') )
    
    if( usr_par.config.n_basis_fct == 0 )
        m_true = zeros( usr_par.config.nx, usr_par.config.nz, 2 );
    else
        m_true = zeros( usr_par.config.nx, usr_par.config.nz, usr_par.config.n_basis_fct+1 );
    end
    
    m_true(:,:,1:end-1) = make_noise_source( 'gaussian', usr_par.config.n_basis_fct );
    m_true(:,:,end) = define_material_parameters( usr_par.config.nx, usr_par.config.nz, 999 );
    
    m_true = map_m_to_parameters( map_parameters_to_m( m_true, usr_par ) , usr_par );
    m_parameters = map_m_to_parameters(model.m, usr_par) - m_true;


%% plot inversion model
else
   
    if( exist('path2', 'var') )
        model2 = load(path2);
        m_parameters = map_m_to_parameters(model.m, usr_par) - map_m_to_parameters(model2.model.m, usr_par);
    else
        m_parameters = map_m_to_parameters(model.m, usr_par);
    end
    
end


%% load array
load ../output/interferometry/array_16_ref.mat
if( ~exist('array', 'var') )
    array = [];
end


%% actual plotting commadn
if( exist('clim','var') )
    plot_models( m_parameters, array, [clim(1) clim(2) clim(3) clim(4)] );
else
%     cm = cbrewer('div','RdBu',120,'PCHIP');
%     plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 0 0 0], 'no', 'no', cm );

    load clim.mat
%     plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 0 clim(1) clim(2)], 'no', 'no' );
    plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 4 4.6e10 5.0e10], 'no', 'no' );
%     plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 0 0 0], 'no', 'no' );

end

if( strcmp(gradient, 'yes') )
    
    gradparameters = map_gradm_to_gradparameters( 0, model.gradient, usr_par );
    max_source = max(max( abs( gradparameters(:,:,end-1) ) ));
    max_structure = max(max( abs( gradparameters(:,:,end) ) ));
    
    cm = cbrewer('div','RdBu',120,'PCHIP');
    plot_models( -gradparameters, usr_par.config.n_basis_fct, array, [-max_source max_source -max_structure max_structure], 'no', 'no', cm, cm );
    
end




return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% plot correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = load('~/Desktop/runs/2016_start/data_iugg/data_16_ref_0_1h1g_iugg_newpara.mat');
% initial = load('~/Desktop/runs/inversion_basis_fct/data/data_16_ref_0_uniform1_homogeneous.mat');
t = data.t;

figure
plot_recordings(data.c_data,t,'vel','k',true);
% plot_recordings(initial.c_data,t,'vel','r',true);
plot_recordings(model.correlation,t,'vel','b',true);



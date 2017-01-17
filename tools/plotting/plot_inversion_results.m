
clear all
% close all
% clc

[Lx, Lz, nx, nz, dt, nt, order, model_type, source_type, n_basis_fct, fw_nth] = input_parameters();


% path = '~/Desktop/model_00.mat';
% path = '~/Diss/Paper/phd_paper_1/matlab_files/final_models/joint_model_5.mat';
% path = '~/Diss/Paper/phd_paper_1/matlab_files/final_models/loga_cc_model_27.mat';
% path = 'initial_models/structure_random_0.10_cc_equal_homogeneous_regu_1em2_smooth_5e4.mat';
path = 'initial_models/structure_random_0.07_cc_equal_homogeneous_regu_1em2_smooth_5e4.mat';
% path2 = '~/Desktop/model_137.mat';


difftrue = 'no';
gradient = 'no';



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
    
    
    if( isfield(model, 'sigma') )
        usr_par.kernel.sigma.source = model.sigma.source;
        usr_par.kernel.sigma.structure = model.sigma.structure;
    end
    
    if( isfield(model, 'imfilter') )
        if( isfield(model.imfilter, 'source') )
            usr_par.kernel.imfilter.source = model.imfilter.source;
            usr_par.kernel.imfilter.structure = model.imfilter.structure;
        else
            usr_par.kernel.imfilter.source = model.imfilter;
            usr_par.kernel.imfilter.structure = model.imfilter;
        end
    end
    
    if( isfield(model, 'config') )
        if( isfield(model.config, 'n_basis_fct') )
            usr_par.config.n_basis_fct = model.config.n_basis_fct;
        else
            usr_par.config.n_basis_fct = n_basis_fct;
        end
    else
        usr_par.config.n_basis_fct = n_basis_fct;
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


if( usr_par.config.n_basis_fct > 5 )
    
    [~,n_sample] = input_interferometry();
    dist_general = zeros( nx, nz, n_sample );
    grad_general = zeros( nx, nz, n_sample );
    int_limits = integration_limits( n_sample, usr_par.config.n_basis_fct );
    
    model.gradient = reshape(model.gradient,nx,nz,usr_par.config.n_basis_fct+1);
    
    for is = 1:n_sample
        
        for ib = 1:usr_par.config.n_basis_fct
            if( is >= int_limits(ib,1) && is <= int_limits(ib,2) )
                break;
            end
        end
        
        dist_general(:,:,is) = m_parameters(:,:,ib);
        grad_general(:,:,is) = model.gradient(:,:,ib);
        
    end
    
    dist_few = zeros( nx, nz, 5 );
    grad_few = zeros( nx, nz, 5 );
    int_limits = integration_limits( n_sample, 5 );
    
    for ib = 1:5
        
        indices = int_limits(ib,1) : int_limits(ib,2);
        ni = length(indices);
        
        for k = indices
            dist_few(:,:,ib) = dist_few(:,:,ib) + dist_general(:,:,k) / ni;
            grad_few(:,:,ib) = grad_few(:,:,ib) + grad_general(:,:,k) / ni;
        end
        
    end
    
    tmp = zeros( nx, nz, 6 );
    tmp(:,:,1:5) = dist_few;
    tmp(:,:,6) = m_parameters(:,:,end);
    m_parameters = tmp;
    
    tmp(:,:,1:5) = grad_few;
    tmp(:,:,6) = model.gradient(:,:,end);
    model.gradient = tmp;
    
    usr_par.config.n_basis_fct = 5;
end


%% load array
load ../output/interferometry/array_16_ref.mat
if( ~exist('array', 'var') )
    array = [];
end


%% convert to velocity
% material parameters
if( isempty( find( m_parameters(:,:,end) < 0, 1 )) )
    [~,rho] = define_material_parameters( nx, nz, model_type );
    m_parameters(:,:,end) = sqrt( m_parameters(:,:,end)./rho );
end


%% actual plotting command
if( exist('clim','var') )
    plot_models( m_parameters, array, [clim(1) clim(2) clim(3) clim(4)] );
else
    
    if( ~exist('path2', 'var') )
        
        cm_source_orig = cbrewer('div','RdBu',120,'PCHIP');

        % cm_source = cm_source_orig(13:120,:);
        % plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 2.2 3700 4300], 'no', 'no', cm_source );
        % plot_models_poster( m_parameters, usr_par.config.n_basis_fct, array, [0 2.2 3700 4300], 'no', 'no', cm_source );
        
        cm_source = cm_source_orig(48:120,:);
        % plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 5.3 3700 4300], 'no', 'no', cm_source );
        plot_models_poster( m_parameters, usr_par.config.n_basis_fct, array, [0 5.3 3700 4300], 'no', 'no', cm_source );
        
        % cm_source = cm_source_orig(50:120,:);
        % plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 7 3700 4300], 'no', 'no', cm_source );
        
        % plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 0 3700 4300], 'no', 'no' );
        
    else
        plot_models( m_parameters, usr_par.config.n_basis_fct, array, [0 0 0 0], 'no', 'no' );
    end

end

if( strcmp(gradient, 'yes') )
    
%     grad_parameters = map_gradm_to_gradparameters( 0, model.gradient, usr_par );
    grad_parameters = reshape( model.gradient, nx, nz, [] );
    max_source = max(max( abs( grad_parameters(:,:,end-1) ) ));
    max_structure = max(max( abs( grad_parameters(:,:,end) ) ));
    
    cm = cbrewer('div','RdBu',120,'PCHIP');
    plot_models( -grad_parameters, usr_par.config.n_basis_fct, array, [-max_source max_source -max_structure max_structure], 'no', 'no', cm, cm );
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% plot correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return
% fprintf('\nResidual nach Inversion:\n')
% data = load('../output/interferometry/data_16_ref_0_gaussian_random_0.07_0.8e10_nosmooth.mat');
% misfit_1 = compare_influences( model.correlation, data.c_data, data.t, array, ref_stat, 'no' );

% return
% fprintf('\nEffekt - Struktur:\n')
% data2 = ;
% misfit_2 = compare_influences( data2.c_data, data.c_data, data.t, array, ref_stat, 'no' );

% fprintf('\nEffekt - Struktur und Quelle:\n')
% data3 = ;
% misfit_3 = compare_influences( data3.c_data, data.c_data, data.t, array, ref_stat, 'no' );

return
initial = load('~/Desktop/runs/paper/data_inversion/data_16_ref_0_equal_homogeneous.mat');

figure
data = load('../output/interferometry/data_16_ref_0_gaussian_random_0.07_0.8e10_nosmooth.mat');
plot_recordings( data.c_data, data.t, 'vel', 'k', true );
plot_recordings( initial.c_data, data.t, 'vel', 'r', true );
plot_recordings( model.correlation, data.t, 'vel', 'b', true );



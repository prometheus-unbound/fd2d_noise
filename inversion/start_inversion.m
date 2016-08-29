
addpath(genpath('../'))
[~,~,nx,nz,~,~,~,model_type,source_type,n_basis_fct] = input_parameters();
usr_par.config.nx = nx; 
usr_par.config.nz = nz;
usr_par.config.n_basis_fct = n_basis_fct; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.cluster = 'monch';
% 'local';
% 'monch';
% 'euler';
% 'brutus';


usr_par.type = 'source';
% 'source'
% 'structure'
% 'joint'


usr_par.use_mex = 'yes';
% 'yes'
% 'no'


% define reference models for perturbations and mu_0 for structure
usr_par.initial.ref_source = 0;
usr_par.initial.ref_structure = 1;
usr_par.initial.mu_0 = 4.8e10;


% if wanted, load initial perturbations for source and/or structure
% initial_model = load('initial_models/model_56.mat');
% usr_par.initial.source_dist = initial_model.model.m( 1 : end - nx*nz, 1 );
% initial_model = load('initial_models/structure_cc_random_0.1_norm.mat');
% usr_par.initial.structure = initial_model.model.m( end - nx*nz + 1 : end, 1 );


usr_par.measurement.source = 'amplitude_difference';
usr_par.measurement.structure = 'waveform_difference';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';


% define weighting for kernels - only necessary for joint inversion
usr_par.kernel.weighting = 0.5;


% design gaussian filter for smoothing of kernel (set second input variable to [1,1] to turn it off)
% usr_par.kernel.imfilter.source = fspecial('gaussian',[75 75], 30);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[1 1], 1);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[40 40], 20);
% usr_par.kernel.imfilter.source = fspecial('gaussian',[1 1], 1);
% usr_par.kernel.imfilter.structure = usr_par.kernel.imfilter.source;

usr_par.kernel.sigma.source = [1 1];
usr_par.kernel.sigma.structure = usr_par.kernel.sigma.source;



% parameterize source distribution as ring
usr_par.ring.switch = 'no';
usr_par.ring.x_center_ring = 1.0e6;
usr_par.ring.z_center_ring = 1.0e6;
usr_par.ring.radius = 6.4e5;
usr_par.ring.thickness = 2.0e5;
usr_par.ring.taper_strength = 70e8;


% load array with reference stations and data
usr_par.network = load('../output/interferometry/array_16_ref.mat');
usr_par.data = load('../output/interferometry/data_16_ref_0_h1g_iugg_smooth.mat');


% do measurement on displacement or velocity correlations (for NOW: use 'dis')
usr_par.veldis = 'dis';


% filter correlations for inversion; if yes, specify f_min and f_max
usr_par.filter.apply_filter = 'no';
usr_par.filter.f_min = 1/7 - 0.01;
usr_par.filter.f_max = 1/7 + 0.01;


% regularization j_total = dj/dm + alpha * ||m_source-m0_source||^2 + beta * ||m_structure-m0_structure||^2 
% (i.e. set alpha or beta to zero to turn it off)
usr_par.regularization.alpha = 1e-6;
usr_par.regularization.beta = 0;
usr_par.regularization.weighting = weighting( nx, nz );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FOR LBFGS
options.tolerance = 1e-3;
options.successive_change = 1e-16;
options.successive_iterations = 400;
options.max_iterations = 500;
options.wolfe_try_to_increase_step_length = false;
options.verbose = true;

if( ~strcmp( usr_par.type, 'source' ) && strcmp( usr_par.measurement.structure, 'cc_time_shift') )
    options.init_step_length = 0.025;
elseif( strcmp( usr_par.type, 'source' ) && strcmp( usr_par.measurement.source, 'log_amplitude_ratio') )
    options.init_step_length = 8.0;
else
    options.init_step_length = 1.0;
end


%%% FOR STEEPEST DESCENT
% options.init_step_length = 1.0;
% options.stepsize = 1.0;
% options.verbose = false;
% options.output_file = 'iterations_lbfgs.tab';
% options.max_iterations = 100;
% options.tolerance = 1e-3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% running the inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( strcmp( usr_par.type, 'source') )
    
    if( usr_par.kernel.weighting ~= 0.0 )
        
        usr_par.kernel.weighting = 0.0;
        usr_par.regularization.beta = 0.0;
        
        fprintf('\nset usr_par.kernel.weighting to 0.0')
        fprintf('\nset usr_par.regularization.beta to 0.0\n')
        
    end
    
elseif( strcmp( usr_par.type, 'structure') )
    
    if( usr_par.kernel.weighting ~= 1.0 )
        
        usr_par.kernel.weighting = 1.0;
        usr_par.regularization.alpha = 0.0;
    
        fprintf('\nset usr_par.kernel.weighting to 1.0')
        fprintf('\nset usr_par.regularization.alpha to 0.0\n')
    end
    
end


% start matlabpool
if( ~strcmp( usr_par.cluster, 'local') )
    parobj = start_cluster( usr_par.cluster, '', size(usr_par.network.ref_stat,1));
end


% set necessary fields that might not have been set
[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% set up initial model
if( usr_par.config.n_basis_fct == 0 )
    m_parameters = zeros(nx, nz, 2);
else
    m_parameters = zeros(nx, nz, usr_par.config.n_basis_fct+1);
end

m_parameters(:,:,1:end-1) = make_noise_source( source_type, usr_par.config.n_basis_fct );
m_parameters(:,:,end) = define_material_parameters( nx, nz, model_type );


% convert to optimization variable
m0 = map_parameters_to_m( m_parameters, usr_par );


% for a restart of the inversion m0 is given by model_type and source_type in input_parameters.m
usr_par.m0 = m0;


if( isfield( usr_par, 'initial') )
    if( isfield( usr_par.initial, 'source_dist') )
        m0( 1 : end - nx*nz, 1 ) = usr_par.initial.source_dist;
    end
    
    if( isfield( usr_par.initial, 'structure') )
        m0( end - nx*nz + 1 : end, 1 ) = usr_par.initial.structure;
    end
end


% uncomment next line if the specified initial model is the real initial model and not just a restart of the inversion
% usr_par.m0 = m0;


% run inversion
[flag, mfinal, usr_par] = optlib_lbfgs( m0, options, usr_par );


% save solution
save( '../output/solution.mat', 'flag', 'mfinal', 'usr_par', '-v7.3' )


% close matlabpool
if( ~strcmp( usr_par.cluster, 'local' ) )
    delete(parobj)
end


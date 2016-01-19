
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usr_par.cluster = 'monch';
% 'local';
% 'monch';
% 'euler';
% 'brutus';


usr_par.use_mex = 'yes';


usr_par.type = 'source';
% 'source';
% 'source_constrained';
% 'structure';

if( strcmp( usr_par.type, 'structure') )
    addpath(genpath('../'))
    [~,~,nx,nz,~,~,~,model_type] = input_parameters();
    [~,rho] = define_material_parameters( nx, nz, model_type);

    usr_par.structure_inversion.rho = rho;
    usr_par.structure_inversion.v0 = 4000;    
end


usr_par.measurement = 'log_amplitude_ratio';
% 'log_amplitude_ratio';
% 'amplitude_difference';
% 'waveform_difference';
% 'cc_time_shift';


% do measurement on displacement or velocity correlations (for NOW: use 'dis')
usr_par.veldis = 'dis';


% filter correlations for inversion; if yes, specify f_min and f_max
usr_par.filter.apply_filter = 'no';
usr_par.filter.f_min = 1/7 - 0.01;
usr_par.filter.f_max = 1/7 + 0.01;


% load array with reference stations and data
usr_par.network = load('../output/interferometry/array_16_ref.mat');
usr_par.data = load('../output/interferometry/data_16_ref_0_1h1g_iugg.mat');


% load source distribution that should be used for structure inversion
% usr_par.source_dist = load('initial_models/model_68.mat');


% specify percentile for clipping of kernel      ( = 0 to turn it off )
usr_par.kernel.percentile = 0;


% desing gaussian filter for smoothing of kernel
usr_par.kernel.imfilter = fspecial('gaussian',[50 50], 20);


% debug mode
usr_par.debug.switch = 'no';
usr_par.debug.df = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FOR LBFGS
options.tolerance = 1e-3;
options.successive_change = 1e-16;
options.successive_iterations = 100;
options.max_iterations = 100;
options.wolfe_try_to_increase_step_length = false;
options.verbose = true;
if( strcmp( usr_par.type, 'structure') )
    options.init_step_length = 0.5;
else
    options.init_step_length = 16.0;
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


% start matlabpool and set up path
if( ~strcmp( usr_par.cluster, 'local') )
    addpath(genpath('../'))
    parobj = start_cluster( usr_par.cluster, '', size(usr_par.network.ref_stat,1));
end


% set necessary fields that might not have been set
[usr_par] = usr_par_init_default_parameters_lbfgs(usr_par);


% run source inversion
if( strcmp( usr_par.type, 'source' ) )
    
    % set up initial model
    m0 = make_noise_source();
    m0 = reshape(m0,[],1);
    
    % do inversion
    [flag, mfinal, usr_par] = optlib_lbfgs(m0, options, usr_par);
    
    % save solution
    save('../output/solution.mat', 'flag', 'mfinal', 'usr_par')
    
    
% run source inversion with lower and upper bounds
elseif( strcmp( usr_par.type, 'source_constrained' ) )
    
    usr_par.type = 'source';
    
    % set up initial model and lower/upper bounds
    m0 = make_noise_source();
    m0 = reshape(m0,[],1);    
    xl = 0 * m0;
    xu = inf * m0;
    
    % do inversion
    [m, c_final] = projected_steepest_descent( m0, xl, xu, 'get_obj_grad', 0.05, 0);

    % save solution
    save('../output/solution.mat', 'm', 'c_final')
        
    
% run structure inversion
elseif( strcmp( usr_par.type, 'structure' ) )
    
    % set up initial model
    [~,~,nx,nz] = input_parameters();
    m0 = zeros(nx*nz, 1);
    
    % do inversion
    % [flag, mfinal, usr_par] = optlib_steepest_descent(m0, options, usr_par);
    [flag, mfinal, usr_par] = optlib_lbfgs(m0, options, usr_par);
    
    % save solution
    save('../output/solution.mat', 'flag', 'mfinal', 'usr_par')    
    
end


% close matlabpool and clean up path
if( ~strcmp( usr_par.cluster, 'local' ) )
    delete(parobj)
    rmpath(genpath('../'))
end

